%% === AGRÉGATION STATS (SANS CALCUL) — copie des métriques + 'mean' jambe ===
% Lit chaque Coherence_<ID>.mat dans \ALL, copie:
%   - MeanCoherence_*_*_* (full et par sous-phase)
%   - CoherenceArea_*_*_* (full) et CoherenceArea_<Phase>_*_*_* (par sous-phase)
% dans STATS au format :
%   Condition > Sous-phase(Full/LoadingResponse/MidStance/PreSwing/Swing)
%   > Pair(m1_m2) > Jambe(left/right/mean) > Bande(Alpha/Beta/Gamma)
%   > Metric(mean_coherence, coherence_area)
% Seule la jambe 'mean' est calculée ici (moyenne des côtés dispo).
clear; clc; close all;

% Dossiers
coh_root = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(coh_root, 'ALL');
out_dir  = fullfile(coh_root, 'STATISTIQUE');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% Constantes
bandNames  = {'Alpha','Beta','Gamma'};
phaseNames = {'Full','LoadingResponse','MidStance','PreSwing','Swing'};
isPhase    = @(p) any(strcmp(p, phaseNames(2:end))); % sans 'Full'
toPairKey  = @(m1,m2) sprintf('%s_%s', m1, m2);

% Fichiers
S = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(S)
    error('Aucun fichier "Coherence_*.mat" dans %s', all_dir);
end
fprintf('>> %d fichier(s) trouvés dans \\ALL.\n', numel(S));

for iF = 1:numel(S)
    fpath = fullfile(S(iF).folder, S(iF).name);
    [~, base, ~] = fileparts(fpath);
    tok = regexp(base, '^Coherence_(.+)$', 'tokens', 'once');
    if isempty(tok), warning('Nom inattendu: %s (skip)', base); continue; end
    pid = tok{1};
    fprintf('\n==============================\nParticipant : %s\n', pid);

    L = load(fpath);
    if ~isfield(L, 'DATA'), warning('Pas de DATA dans %s (skip)', fpath); continue; end
    DATA = L.DATA;

    % Détection conditions exploitables
    condNames = fieldnames(DATA);
    keep = false(size(condNames));
    for k = 1:numel(condNames)
        nm = condNames{k};
        keep(k) = isstruct(DATA.(nm)) && (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
    end
    condNames = condNames(keep);
    if isempty(condNames), warning('Aucune condition exploitable pour %s', pid); continue; end

    STATS = struct();
    STATS.meta.participant_id = pid;
    STATS.meta.metrics = {'mean_coherence','coherence_area'};
    STATS.meta.bands   = bandNames;
    STATS.meta.phases  = phaseNames;
    notes = {};

    for iC = 1:numel(condNames)
        cond = condNames{iC};
        sides_present = intersect({'left','right'}, fieldnames(DATA.(cond))');

        % Cache pour construire la jambe "mean" ensuite
        cache = struct(); % cache.(phase).(pair).(band).left/right.(metric)

        for s = 1:numel(sides_present)
            sd = sides_present{s};
            Sside = DATA.(cond).(sd);
            fns   = fieldnames(Sside);

            % On s'appuie UNIQUEMENT sur les champs MeanCoherence_* pour répertorier
            for iFN = 1:numel(fns)
                fn = fns{iFN};
                if ~startsWith(fn,'MeanCoherence_'), continue; end

                % Cas 1 : full phase -> MeanCoherence_<Band>_<m1>_<m2>
                Tfull = regexp(fn, '^MeanCoherence_(?<band>Alpha|Beta|Gamma)_(?<m1>[^_]+)_(?<m2>[^_]+)$', 'names');
                % Cas 2 : phase -> MeanCoherence_<Phase>_<Band>_<m1>_<m2>
                Tph   = regexp(fn, '^MeanCoherence_(?<phase>LoadingResponse|MidStance|PreSwing|Swing)_(?<band>Alpha|Beta|Gamma)_(?<m1>[^_]+)_(?<m2>[^_]+)$', 'names');

                if ~isempty(Tfull)
                    phaseKey = 'Full';
                    band  = Tfull.band; m1 = Tfull.m1; m2 = Tfull.m2;
                elseif ~isempty(Tph)
                    phaseKey = Tph.phase;
                    band  = Tph.band;   m1 = Tph.m1;   m2 = Tph.m2;
                else
                    continue; % autre format non concerné
                end

                pairKey = toPairKey(m1,m2);
                mc_val  = Sside.(fn); % déjà calculé par le 1er script

                % Récupérer l'aire correspondante (sans recalcul)
                if strcmp(phaseKey,'Full')
                    area_field = sprintf('CoherenceArea_%s_%s_%s', band, m1, m2);
                else
                    area_field = sprintf('CoherenceArea_%s_%s_%s_%s', phaseKey, band, m1, m2);
                end
                if isfield(Sside, area_field)
                    area_val = Sside.(area_field);
                else
                    area_val = NaN;
                    notes{end+1} = sprintf('%s | %s | %s: %s manquant', cond, sd, pairKey, area_field);
                end

                % Écriture (copie) dans STATS (aucun calcul)
                STATS.(cond).(phaseKey).(pairKey).(sd).(band).mean_coherence = mc_val;
                STATS.(cond).(phaseKey).(pairKey).(sd).(band).coherence_area = area_val;

                % Mise en cache pour 'mean' jambe
                cache.(phaseKey).(pairKey).(band).(sd).mean_coherence = mc_val;
                cache.(phaseKey).(pairKey).(band).(sd).coherence_area = area_val;
            end
        end

        % Construire la jambe "mean" (moyenne des côtés disponibles, NaN-robuste)
        if ~isempty(fieldnames(cache))
            phases = fieldnames(cache);
            for ip = 1:numel(phases)
                ph = phases{ip};
                pairs = fieldnames(cache.(ph));
                for ipa = 1:numel(pairs)
                    pk = pairs{ipa};
                    for ib = 1:numel(bandNames)
                        b = bandNames{ib};

                        vals_mc = [];
                        vals_ar = [];

                        if isfield(cache.(ph).(pk).(b),'left')
                            vals_mc(end+1) = cache.(ph).(pk).(b).left.mean_coherence; %#ok<AGROW>
                            vals_ar(end+1) = cache.(ph).(pk).(b).left.coherence_area; %#ok<AGROW>
                        end
                        if isfield(cache.(ph).(pk).(b),'right')
                            vals_mc(end+1) = cache.(ph).(pk).(b).right.mean_coherence; %#ok<AGROW>
                            vals_ar(end+1) = cache.(ph).(pk).(b).right.coherence_area; %#ok<AGROW>
                        end

                        if ~isempty(vals_mc) || ~isempty(vals_ar)
                            STATS.(cond).(ph).(pk).mean.(b).mean_coherence = mean(vals_mc, 'omitnan');
                            STATS.(cond).(ph).(pk).mean.(b).coherence_area = mean(vals_ar, 'omitnan');
                        end
                    end
                end
            end
        end
    end

    if ~isempty(notes)
        STATS.meta.notes = unique(notes);
        fprintf('  Notes consignées: %d\n', numel(STATS.meta.notes));
    end

    % Sauvegarde
    out_file = fullfile(out_dir, ['STATS_' pid '.mat']);
    save(out_file, 'STATS', '-v7.3');
    fprintf('  -> Sauvé : %s\n', out_file);
end

fprintf('\n>> Terminé.\n');