%% AGRÉGATION DES MÉTRIQUES DE COHÉRENCE INDIVIDUELLE DANS UNE STRUCTURE STATS
%
% OBJECTIF :
%   Ce script lit les fichiers de résultats de cohérence (Coherence_<PID>.mat,
%   produits par les scripts 7 et 9) et en extrait, pour chaque participant,
%   les métriques déjà calculées sans les recalculer :
%     - Cohérence moyenne par bande (MeanCoherence_*)
%     - Aire sous la courbe significative (CoherenceArea_*)
%   Ces valeurs sont réorganisées dans une structure STATS hiérarchique,
%   indexée par condition → sous-phase → paire → jambe → bande → métrique.
%
%   Une jambe virtuelle 'mean' est calculée ici comme la moyenne des côtés
%   gauche et droit disponibles (robuste aux NaN).
%
%   Un fichier STATS_<PID>.mat est produit par participant dans le dossier
%   STATISTIQUE, prêt pour les analyses de groupe (script 11).
%
% ENTRÉES :
%   - Coherence_<PID>.mat : Fichiers issus des scripts 7 + 9, contenant la
%                           structure DATA avec les champs :
%       * MeanCoherence_<Bande>_<m1>_<m2>                    (cycle complet)
%       * MeanCoherence_<Phase>_<Bande>_<m1>_<m2>            (par sous-phase)
%       * CoherenceArea_<Bande>_<m1>_<m2>                    (cycle complet)
%       * CoherenceArea_<Phase>_<Bande>_<m1>_<m2>            (par sous-phase)
%
% SORTIES :
%   - STATS_<PID>.mat : Un fichier par participant, contenant la structure
%                       STATS organisée comme suit :
%       STATS.(cond).(phase).(pair).(side).(band).mean_coherence
%       STATS.(cond).(phase).(pair).(side).(band).coherence_area
%       où side ∈ {left, right, mean}
%            phase ∈ {Full, LoadingResponse, MidStance, PreSwing, Swing}
%            band  ∈ {Alpha, Beta, Gamma}
%       * STATS.meta : métadonnées (participant, bandes, phases, notes)
%
% NOTES :
%   - Aucun recalcul de cohérence n'est effectué ici : les valeurs sont
%     simplement copiées depuis DATA et réorganisées.
%   - La jambe 'mean' est la seule valeur calculée (moyenne NaN-robuste).
%   - Les champs CoherenceArea absents génèrent un NaN et une note dans
%     STATS.meta.notes.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ===================== CHEMINS & CONSTANTES =============================

% Dossiers d'entrée et de sortie
coh_root = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(coh_root, 'ALL');       % Fichiers Coherence_<PID>.mat
out_dir  = fullfile(coh_root, 'STATISTIQUE'); % Fichiers STATS_<PID>.mat de sortie
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

% Constantes — doivent correspondre exactement aux suffixes des champs DATA
bandNames  = {'Alpha','Beta','Gamma'};
phaseNames = {'Full','LoadingResponse','MidStance','PreSwing','Swing'};

% Fonction anonyme : construit la clé d'une paire 'muscle1_muscle2'
toPairKey  = @(m1,m2) sprintf('%s_%s', m1, m2);

%% ===================== DÉTECTION DES FICHIERS ===========================

S = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(S)
    error('Aucun fichier "Coherence_*.mat" dans %s', all_dir);
end
fprintf('>> %d fichier(s) trouvés dans \\ALL.\n', numel(S));

%% ===================== BOUCLE PRINCIPALE PAR PARTICIPANT ================

for iF = 1:numel(S)
    fpath = fullfile(S(iF).folder, S(iF).name);

    % Extraction de l'identifiant du participant depuis le nom de fichier
    % (ex: 'Coherence_CTL_63.mat' → 'CTL_63')
    [~, base, ~] = fileparts(fpath);
    tok = regexp(base, '^Coherence_(.+)$', 'tokens', 'once');
    if isempty(tok)
        warning('Nom de fichier inattendu: %s (ignoré)', base);
        continue;
    end
    pid = tok{1};
    fprintf('\n==============================\nParticipant : %s\n', pid);

    % Chargement de la variable DATA
    L = load(fpath);
    if ~isfield(L, 'DATA')
        warning('Pas de variable DATA dans %s (ignoré)', fpath);
        continue;
    end
    DATA = L.DATA;

    % --- Filtrage des conditions exploitables ---
    % Conserve uniquement les champs de DATA qui contiennent au moins
    % un côté ('left' ou 'right') avec des résultats de cohérence
    condNames = fieldnames(DATA);
    keep = false(size(condNames));
    for k = 1:numel(condNames)
        nm = condNames{k};
        keep(k) = isstruct(DATA.(nm)) && ...
                  (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
    end
    condNames = condNames(keep);
    if isempty(condNames)
        warning('Aucune condition exploitable pour %s', pid);
        continue;
    end

    % --- Initialisation de la structure STATS du participant ---
    STATS = struct();
    STATS.meta.participant_id = pid;
    STATS.meta.metrics        = {'mean_coherence','coherence_area'};
    STATS.meta.bands          = bandNames;
    STATS.meta.phases         = phaseNames;
    notes = {}; % Journal des avertissements non bloquants

    for iC = 1:numel(condNames)
        cond = condNames{iC};

        % Côtés effectivement présents pour cette condition
        sides_present = intersect({'left','right'}, fieldnames(DATA.(cond))');

        % Cache intermédiaire pour construire la jambe 'mean' après la boucle
        % Structure : cache.(phase).(pair).(band).(side).(metric)
        cache = struct();

        for s = 1:numel(sides_present)
            sd    = sides_present{s};
            Sside = DATA.(cond).(sd);
            fns   = fieldnames(Sside);

            % --- Parcours des champs MeanCoherence_* ---
            % On utilise uniquement ces champs comme pivot pour identifier
            % les paires et sous-phases disponibles (évite de chercher
            % directement les champs CoherenceArea qui pourraient manquer)
            for iFN = 1:numel(fns)
                fn = fns{iFN};
                if ~startsWith(fn,'MeanCoherence_'), continue; end

                % Tentative de correspondance avec les deux formats possibles :
                %   Format 1 (cycle complet) : MeanCoherence_<Band>_<m1>_<m2>
                %   Format 2 (sous-phase)    : MeanCoherence_<Phase>_<Band>_<m1>_<m2>
                Tfull = regexp(fn, ...
                    '^MeanCoherence_(?<band>Alpha|Beta|Gamma)_(?<m1>[^_]+)_(?<m2>[^_]+)$', ...
                    'names');
                Tph = regexp(fn, ...
                    '^MeanCoherence_(?<phase>LoadingResponse|MidStance|PreSwing|Swing)_(?<band>Alpha|Beta|Gamma)_(?<m1>[^_]+)_(?<m2>[^_]+)$', ...
                    'names');

                if ~isempty(Tfull)
                    phaseKey = 'Full';
                    band = Tfull.band; m1 = Tfull.m1; m2 = Tfull.m2;
                elseif ~isempty(Tph)
                    phaseKey = Tph.phase;
                    band = Tph.band;   m1 = Tph.m1;   m2 = Tph.m2;
                else
                    continue; % Champ MeanCoherence non reconnu → ignoré
                end

                pairKey = toPairKey(m1, m2);
                mc_val  = Sside.(fn); % Valeur de cohérence moyenne (déjà calculée)

                % Récupération de l'aire correspondante (champ produit par script 9)
                % Si absent : NaN + note dans le journal
                if strcmp(phaseKey,'Full')
                    area_field = sprintf('CoherenceArea_%s_%s_%s', band, m1, m2);
                else
                    area_field = sprintf('CoherenceArea_%s_%s_%s_%s', phaseKey, band, m1, m2);
                end

                if isfield(Sside, area_field)
                    area_val = Sside.(area_field);
                else
                    area_val = NaN;
                    notes{end+1} = sprintf('%s | %s | %s: %s manquant', ...
                        cond, sd, pairKey, area_field);
                end

                % Copie des métriques dans STATS (aucun calcul supplémentaire)
                STATS.(cond).(phaseKey).(pairKey).(sd).(band).mean_coherence = mc_val;
                STATS.(cond).(phaseKey).(pairKey).(sd).(band).coherence_area = area_val;

                % Mise en cache pour le calcul ultérieur de la jambe 'mean'
                cache.(phaseKey).(pairKey).(band).(sd).mean_coherence = mc_val;
                cache.(phaseKey).(pairKey).(band).(sd).coherence_area = area_val;
            end
        end

        % --- Calcul de la jambe virtuelle 'mean' ---
        % Moyenne NaN-robuste des côtés disponibles (gauche et/ou droit)
        % pour chaque combinaison phase × paire × bande
        if ~isempty(fieldnames(cache))
            phases = fieldnames(cache);
            for ip = 1:numel(phases)
                ph = phases{ip};
                pairs = fieldnames(cache.(ph));
                for ipa = 1:numel(pairs)
                    pk = pairs{ipa};
                    for ib = 1:numel(bandNames)
                        b = bandNames{ib};

                        % Collecte des valeurs disponibles selon les côtés présents
                        vals_mc = [];
                        vals_ar = [];
                        if isfield(cache.(ph).(pk).(b), 'left')
                            vals_mc(end+1) = cache.(ph).(pk).(b).left.mean_coherence; 
                            vals_ar(end+1) = cache.(ph).(pk).(b).left.coherence_area; 
                        end
                        if isfield(cache.(ph).(pk).(b), 'right')
                            vals_mc(end+1) = cache.(ph).(pk).(b).right.mean_coherence; 
                            vals_ar(end+1) = cache.(ph).(pk).(b).right.coherence_area; 
                        end

                        % Enregistrement de la moyenne bilatérale
                        if ~isempty(vals_mc) || ~isempty(vals_ar)
                            STATS.(cond).(ph).(pk).mean.(b).mean_coherence = mean(vals_mc, 'omitnan');
                            STATS.(cond).(ph).(pk).mean.(b).coherence_area = mean(vals_ar, 'omitnan');
                        end
                    end
                end
            end
        end
    end % Conditions

    % --- Enregistrement des notes dans les métadonnées ---
    if ~isempty(notes)
        STATS.meta.notes = unique(notes);
        fprintf('  Notes consignées: %d\n', numel(STATS.meta.notes));
    end

    % --- Sauvegarde du fichier STATS du participant ---
    out_file = fullfile(out_dir, ['STATS_' pid '.mat']);
    save(out_file, 'STATS', '-v7.3');
    fprintf('  -> Sauvé : %s\n', out_file);
end

fprintf('\n>> Terminé.\n');