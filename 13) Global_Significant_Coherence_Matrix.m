%% === ADDITION POINT-À-POINT DES COHÉRENCES SIGNIFICATIVES — PAR GROUPE D'ÂGE ===
% Pour chaque groupe d'âge (Adultes, Adolescents, Enfants, JeunesEnfants):
%   - filtre les fichiers Coherence_<ID>.mat dont l'ID appartient au groupe
%   - agrège point-à-point (Sum, N, Presence_Count) par cond/side/pair
%   - construit le side 'both' = left + right (matrices & compteurs)
%   - sauvegarde un .mat par groupe: Significant_Coherence_SUM_<Groupe>.mat
%   - Construction de Heatmaps en 2nd partie de script
%
% Nécessite un fichier ParticipantGroup.m qui définit une struct Group avec:
%   Group.Adultes = {'CTL_04', ...}; etc.

clear; clc; close all;

% ===== Paramètres =====
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
participant_group_file = 'C:\Users\defsil00\Documents\Script\ParticipantGroup.m';
assert(isfolder(all_dir), 'Dossier introuvable: %s', all_dir);

% Groupes à traiter (laisser vide [] ou {} pour tous les groupes)
groups_to_run = {};  % ex: {'JeunesEnfants' | 'Enfants' | 'Adolescents' | 'Adultes'}

% ===== Charger la définition des groupes =====
assert(exist(participant_group_file,'file')==2, ...
    'Fichier ParticipantGroup.m introuvable: %s', participant_group_file);
run(participant_group_file);  % doit créer la variable Group (struct)
assert(exist('Group','var')==1 && isstruct(Group), ...
    'ParticipantGroup.m doit définir une variable struct "Group".');

all_groups = fieldnames(Group);
if isempty(groups_to_run)
    groups_to_run = all_groups;
else
    miss = setdiff(groups_to_run, all_groups);
    if ~isempty(miss)
        warning('Groupes inconnus ignorés: %s', strjoin(miss, ', '));
        groups_to_run = setdiff(groups_to_run, miss);
    end
end
if isempty(groups_to_run)
    warning('Aucun groupe valide à traiter.'); return;
end

% ===== Lister les fichiers dispo =====
files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(files)
    fprintf('Aucun fichier Coherence_*.mat trouvé dans: %s\n', all_dir); return;
end
fprintf('>> %d fichier(s) détecté(s) dans %s\n', numel(files), all_dir);

% Helpers
get_id = @(nm) regexp(nm, '^Coherence_([^\.]+)\.mat$', 'tokens', 'once');
isFullCycleSignif = @(fn) (startsWith(fn,'CoherenceSignif_') && numel(strsplit(fn,'_'))==3);
sides = {'left','right'};

% ===== Boucle sur les groupes =====
for ig = 1:numel(groups_to_run)
    Gname = groups_to_run{ig};
    wanted_ids = upper(string(Group.(Gname)));

    % --- Sélectionner les fichiers du groupe ---
    idx_files = [];
    file_ids  = strings(0);
    for i = 1:numel(files)
        tok = get_id(files(i).name);
        if isempty(tok), continue; end
        pid = upper(string(tok{1}));
        if any(strcmp(pid, wanted_ids))
            idx_files(end+1) = i;
            file_ids(end+1)  = pid;
        end
    end

    fprintf('\n================= GROUPE: %s =================\n', Gname);
    if isempty(idx_files)
        fprintf('   (aucun participant trouvé dans ce dossier pour ce groupe)\n');
        continue;
    else
        fprintf('   Participants trouvés (%d): %s\n', numel(idx_files), strjoin(cellstr(file_ids), ', '));
        present_ids = unique(file_ids);
        missing_ids = setdiff(wanted_ids, present_ids);
        if ~isempty(missing_ids)
            fprintf('   (absents dans le dossier: %s)\n', strjoin(cellstr(missing_ids), ', '));
        end
    end

    % === Structures d'agrégation POUR CE GROUPE ===
    Significant_Coherence   = struct();
    Significant_Coherence_N = struct();
    Presence_Count          = struct();
    notes = {};

    % ===== Boucle fichiers (du groupe) =====
    for kf = 1:numel(idx_files)
        iFile = idx_files(kf);
        fpath = fullfile(files(iFile).folder, files(iFile).name);
        fprintf('\n  -> [%d/%d] %s\n', kf, numel(idx_files), files(iFile).name);

        try
            S = load(fpath);
            if ~isfield(S, 'DATA')
                notes{end+1} = sprintf('%s : variable DATA absente', files(iFile).name);
                continue;
            end
            DATA = S.DATA;

            condNames = fieldnames(DATA);
            mask = false(size(condNames));
            for k = 1:numel(condNames)
                nm = condNames{k};
                mask(k) = isstruct(DATA.(nm)) && (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
            end
            condNames = condNames(mask);

            for iC = 1:numel(condNames)
                condName = condNames{iC};
                for s = 1:numel(sides)
                    sideStr = sides{s};
                    if ~isfield(DATA.(condName), sideStr), continue; end
                    Str = DATA.(condName).(sideStr);

                    fns = fieldnames(Str);
                    for iF = 1:numel(fns)
                        fn = fns{iF};
                        if ~isFullCycleSignif(fn), continue; end

                        parts = strsplit(fn,'_');              % {'CoherenceSignif','TAprox','TAdist'}
                        pairName = canon_pair(parts{2}, parts{3}); % Ordre des muscles variable

                        M = Str.(fn);
                        if isempty(M) || ndims(M) ~= 2
                            notes{end+1} = sprintf('%s | %s | %s | %s : matrice vide/non-2D', ...
                                files(iFile).name, condName, sideStr, fn);
                            continue;
                        end

                        % === Correctifs NaN / présence ===
                        if ~isfloat(M), M = double(M); end
                        mask_pos = (~isnan(M)) & (M > 0);  % présence réelle >0, sans NaN
                        M(isnan(M)) = 0;                   % empêche NaN de polluer la somme
                        % Option si tu veux empêcher les annulations par valeurs négatives :
                        % M = max(M, 0);

                        % Créer noeuds si besoin
                        if ~isfield(Significant_Coherence, condName)
                            Significant_Coherence.(condName)   = struct();
                            Significant_Coherence_N.(condName) = struct();
                        end
                        if ~isfield(Significant_Coherence.(condName), sideStr)
                            Significant_Coherence.(condName).(sideStr)   = struct();
                            Significant_Coherence_N.(condName).(sideStr) = struct();
                        end
                        if ~isfield(Presence_Count, condName)
                            Presence_Count.(condName) = struct();
                        end
                        if ~isfield(Presence_Count.(condName), sideStr)
                            Presence_Count.(condName).(sideStr) = struct();
                        end

                        % INIT ou ADD
                        if ~isfield(Significant_Coherence.(condName).(sideStr), pairName)
                            Significant_Coherence.(condName).(sideStr).(pairName)   = M;
                            Significant_Coherence_N.(condName).(sideStr).(pairName) = 1;
                            Presence_Count.(condName).(sideStr).(pairName)          = double(mask_pos);

                            % Stocker Freq une fois si dispo (Hz)
                            if ~isfield(Significant_Coherence,'meta') || ~isfield(Significant_Coherence.meta,'Freq')
                                if isfield(Str,'Freq') && ~isempty(Str.Freq) && numel(Str.Freq)==size(M,1)
                                    Significant_Coherence.meta.Freq = Str.Freq(:);
                                end
                            end
                            fprintf('     [+] init %s | %s | %s (%s)\n', condName, sideStr, pairName, mat2str(size(M)));
                        else
                            Msum = Significant_Coherence.(condName).(sideStr).(pairName);
                            if ~isequal(size(Msum), size(M))
                                notes{end+1} = sprintf('%s | %s | %s | %s : taille différente (existant %s vs courant %s) -> SKIP', ...
                                    files(iFile).name, condName, sideStr, pairName, mat2str(size(Msum)), mat2str(size(M)));
                                continue;
                            end
                            Significant_Coherence.(condName).(sideStr).(pairName) = Msum + M;
                            Significant_Coherence_N.(condName).(sideStr).(pairName) = ...
                                Significant_Coherence_N.(condName).(sideStr).(pairName) + 1;

                            % Ajout présence (compte de sujets "positifs" par (f,t))
                            Presence_Count.(condName).(sideStr).(pairName) = ...
                                Presence_Count.(condName).(sideStr).(pairName) + double(mask_pos);
                        end
                    end
                end
            end

        catch ME
            notes{end+1} = sprintf('ERREUR %s : %s', files(iFile).name, ME.message);
            continue;
        end
    end

    % ===== Récap par groupe =====
    fprintf('\n----- RÉCAP %s -----\n', Gname);
    condList = setdiff(fieldnames(Significant_Coherence), {'meta'});
    for iC = 1:numel(condList)
        condName = condList{iC};
        fprintf('Condition: %s\n', condName);
        for s = ["left","right"]
            sideStr = char(s);
            if ~isfield(Significant_Coherence.(condName), sideStr), continue; end
            pairList = fieldnames(Significant_Coherence.(condName).(sideStr));
            fprintf('  %s: %d paire(s)\n', sideStr, numel(pairList));
        end
    end
    if ~isempty(notes)
        fprintf('\nNotes/Warnings (%d):\n', numel(notes));
        max_show = min(40, numel(notes));
        for i = 1:max_show, fprintf('  - %s\n', notes{i}); end
        if numel(notes) > max_show, fprintf('  ... (+%d supplémentaires)\n', numel(notes)-max_show); end
    end

    % ===== Construire 'both' (somme gauche+droite) =====
    condList = setdiff(fieldnames(Significant_Coherence), {'meta'});
    for iC = 1:numel(condList)
        cond = condList{iC};
        hasL = isfield(Significant_Coherence.(cond),'left');
        hasR = isfield(Significant_Coherence.(cond),'right');
        if ~hasL && ~hasR, continue; end

        Lpairs = {}; Rpairs = {};
        if hasL, Lpairs = fieldnames(Significant_Coherence.(cond).left);  end
        if hasR, Rpairs = fieldnames(Significant_Coherence.(cond).right); end
        allPairs = unique([Lpairs; Rpairs]);

        if ~isfield(Significant_Coherence.(cond),'both')
            Significant_Coherence.(cond).both = struct();
        end
        if ~isfield(Significant_Coherence_N, cond), Significant_Coherence_N.(cond) = struct(); end
        if ~isfield(Significant_Coherence_N.(cond),'both')
            Significant_Coherence_N.(cond).both = struct();
        end
        if ~isfield(Presence_Count, cond), Presence_Count.(cond) = struct(); end
        if ~isfield(Presence_Count.(cond),'both')
            Presence_Count.(cond).both = struct();
        end

        for ip = 1:numel(allPairs)
            pair = allPairs{ip};
            hasL_pair = hasL && isfield(Significant_Coherence.(cond).left,  pair);
            hasR_pair = hasR && isfield(Significant_Coherence.(cond).right, pair);
            if ~hasL_pair && ~hasR_pair
    continue;
            end

            if hasL_pair && hasR_pair
                A = Significant_Coherence.(cond).left.(pair);
                B = Significant_Coherence.(cond).right.(pair);
                if ~isequal(size(A), size(B))
                    notes{end+1} = sprintf('%s | both(sum) | %s : tailles L/R diff (%s vs %s) -> SKIP', ...
                        cond, pair, mat2str(size(A)), mat2str(size(B)));
                    continue;
                end
                % Sécuriser contre NaN résiduels (par prudence)
                A(isnan(A)) = 0; B(isnan(B)) = 0;
                Msum = A + B;

                NL = 0; if isfield(Significant_Coherence_N.(cond),'left')  && isfield(Significant_Coherence_N.(cond).left, pair),  NL = Significant_Coherence_N.(cond).left.(pair);  end
                NR = 0; if isfield(Significant_Coherence_N.(cond),'right') && isfield(Significant_Coherence_N.(cond).right,pair), NR = Significant_Coherence_N.(cond).right.(pair); end
                Nsum = NL + NR;

                if isfield(Presence_Count.(cond),'left') && isfield(Presence_Count.(cond).left, pair) && ...
                   isfield(Presence_Count.(cond),'right') && isfield(Presence_Count.(cond).right, pair)
                    Pl = Presence_Count.(cond).left.(pair);
                    Pr = Presence_Count.(cond).right.(pair);
                    Pl(isnan(Pl)) = 0; Pr(isnan(Pr)) = 0;
                    Psum = Pl + Pr;
                else
                    Psum = [];
                end

            elseif hasL_pair
                Msum = Significant_Coherence.(cond).left.(pair);
                Msum(isnan(Msum)) = 0;
                Nsum = 0; if isfield(Significant_Coherence_N.(cond),'left') && isfield(Significant_Coherence_N.(cond).left, pair)
                    Nsum = Significant_Coherence_N.(cond).left.(pair);
                end
                if isfield(Presence_Count.(cond),'left') && isfield(Presence_Count.(cond).left, pair)
                    Psum = Presence_Count.(cond).left.(pair); Psum(isnan(Psum)) = 0;
                else
                    Psum = [];
                end

            else % hasR_pair seulement
                Msum = Significant_Coherence.(cond).right.(pair);
                Msum(isnan(Msum)) = 0;
                Nsum = 0; if isfield(Significant_Coherence_N.(cond),'right') && isfield(Significant_Coherence_N.(cond).right, pair)
                    Nsum = Significant_Coherence_N.(cond).right.(pair);
                end
                if isfield(Presence_Count.(cond),'right') && isfield(Presence_Count.(cond).right, pair)
                    Psum = Presence_Count.(cond).right.(pair); Psum(isnan(Psum)) = 0;
                else
                    Psum = [];
                end
            end

            Significant_Coherence.(cond).both.(pair)   = Msum;
            Significant_Coherence_N.(cond).both.(pair) = Nsum;
            if ~isempty(Psum), Presence_Count.(cond).both.(pair) = Psum; end
        end
    end

    % ===== Sauvegarde par groupe =====
    out_file = fullfile(all_dir, sprintf('Significant_Coherence_SUM_%s.mat', Gname));
    save(out_file, 'Significant_Coherence', 'Significant_Coherence_N', 'Presence_Count', '-v7.3');
    fprintf('>> Sauvegardé: %s\n', out_file);
end

fprintf('\n== Terminé pour groupes: %s ==\n', strjoin(groups_to_run, ', '));

%% === HEATMAPS (Sum & Proportion) PAR GROUPE — side=both ====================
clear; clc; close all;

% ---- Paramètres ----
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
out_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\Visualisation_Coherence';
side    = 'both';
FMAX    = 80;     % fréquence max affichée (Hz)
dpi     = 300;
fmt     = 'png';  % 'png' | 'tiff' | 'jpg'

% (Optionnel) limiter à certains groupes
groups_to_plot = {};  % ex: {'Adultes','Adolescents'}

% ---- Détection des fichiers groupe ----
gfiles = dir(fullfile(all_dir, 'Significant_Coherence_SUM_*.mat'));
gfiles = gfiles(~strcmp({gfiles.name}, 'Significant_Coherence_SUM.mat'));

if ~isempty(groups_to_plot)
    keep = false(1, numel(gfiles));
    for i=1:numel(gfiles)
        tok = regexp(gfiles(i).name, '^Significant_Coherence_SUM_(.+)\.mat$', 'tokens','once');
        if ~isempty(tok)
            keep(i) = any(strcmpi(tok{1}, groups_to_plot));
        end
    end
    gfiles = gfiles(keep);
end

if isempty(gfiles)
    error('Aucun fichier groupe trouvé dans %s (pattern Significant_Coherence_SUM_*.mat).', all_dir);
end

% ---- Dossier racine de sortie ----
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% ============================ BOUCLE GROUPES ===============================
for ig = 1:numel(gfiles)
    tok = regexp(gfiles(ig).name, '^Significant_Coherence_SUM_(.+)\.mat$', 'tokens','once');
    if isempty(tok), fprintf('[skip] Fichier inattendu: %s\n', gfiles(ig).name); continue; end
    Gname = tok{1};
    fprintf('\n=== GROUPE: %s ===\n', Gname);

    S = load(fullfile(gfiles(ig).folder, gfiles(ig).name));
    if ~isfield(S, 'Significant_Coherence') || ~isfield(S, 'Significant_Coherence_N') || ~isfield(S, 'Presence_Count')
        fprintf('  [skip] Structures manquantes dans %s\n', gfiles(ig).name); continue;
    end
    Significant_Coherence   = S.Significant_Coherence;
    Significant_Coherence_N = S.Significant_Coherence_N;
    Presence_Count          = S.Presence_Count;

    % Vecteur de fréquences (Hz) si présent
    hasHz = isfield(Significant_Coherence,'meta') && isfield(Significant_Coherence.meta,'Freq');
    if hasHz
        f_all = Significant_Coherence.meta.Freq(:);
    else
        warning('(%s) meta.Freq manquant : l’axe Y sera en indices (bins).', Gname);
    end

    % Dossier de sortie du groupe
    out_group = fullfile(out_dir, Gname);
    if ~exist(out_group, 'dir'), mkdir(out_group); end

    % ---- Boucle conditions ----
    condList = setdiff(fieldnames(Significant_Coherence), {'meta'});
    for ic = 1:numel(condList)
        cond = condList{ic};
        if ~isfield(Significant_Coherence.(cond), side), continue; end

        pairList = fieldnames(Significant_Coherence.(cond).(side));
        fprintf('  %s : %d paires\n', cond, numel(pairList));

        for ip = 1:numel(pairList)
            pair = pairList{ip};

            if ~isfield(Significant_Coherence_N.(cond), side) || ...
               ~isfield(Significant_Coherence_N.(cond).(side), pair) || ...
               ~isfield(Presence_Count.(cond), side) || ...
               ~isfield(Presence_Count.(cond).(side), pair)
                fprintf('   [skip] Données incomplètes pour %s | %s | %s\n', cond, side, pair);
                continue;
            end

            SumMap = Significant_Coherence.(cond).(side).(pair);   % somme d'excès
            % Protection additionnelle (au cas où)
            SumMap(isnan(SumMap)) = 0;

            N      = Significant_Coherence_N.(cond).(side).(pair); % # contributeurs (scalaire)
            P      = Presence_Count.(cond).(side).(pair);          % # sujets>0 par (f,t)

            [nFreq,nTime] = size(SumMap);
            t = linspace(0,100,nTime);

            % Fréquences croissantes + coupe ≤FMAX
            if hasHz && numel(f_all)==nFreq
                f = f_all(:);
            else
                f = (1:nFreq)';   % fallback: index de bin
            end
            [fc, idx] = sort(f, 'ascend');
            SumMap_c  = SumMap(idx, :);
            P_c       = P(idx, :);

            if hasHz
                mask = fc <= FMAX;
            else
                mask = true(size(fc));
            end
            fvis = fc(mask);
            SumV = SumMap_c(mask, :);
            Pv   = P_c(mask, :);
            Prop = Pv ./ max(N,1);     % 0..1

            % === HEATMAP 1 : Sum of excess ===
            fig1 = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
            % Ultime garde-fou
            SumV(isnan(SumV)) = 0;
            imagesc(t, fvis, SumV); set(gca,'YDir','normal'); axis tight; colorbar;
            xlabel('Gait cycle (%)');
            ylabel( tern(hasHz,'Frequency (Hz)','Frequency (bin)') );
            title(sprintf('%s | %s | %s — Sum of excess — %s', cond, side, strrep(pair,'_','-'), FMAX, Gname));
            colormap(parula);
            if hasHz
                yb = [8 12 13 30 31 60]; yb = yb(yb>=min(fvis) & yb<=max(fvis));
                if ~isempty(yb), hold on; yline(yb,'w--'); end
            end
            file1 = fullfile(out_group, sprintf('Sum_%s_%s_%s_%s_Fmax%g.%s', Gname, cond, side, pair, FMAX, fmt));
            exportgraphics(fig1, file1, 'Resolution', dpi);

            % === HEATMAP 2 : Proportion of subjects ===
            fig2 = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
            imagesc(t, fvis, Prop); set(gca,'YDir','normal'); axis tight; colorbar; caxis([0 1]);
            xlabel('Gait cycle (%)');
            ylabel( tern(hasHz,'Frequency (Hz)','Frequency (bin)') );
            title(sprintf('%s | %s | %s — Proportion — %s', cond, side, strrep(pair,'_','-'), FMAX, Gname));
            colormap(hot);
            if hasHz
                yb = [8 12 13 30 31 60]; yb = yb(yb>=min(fvis) & yb<=max(fvis));
                if ~isempty(yb), hold on; yline(yb,'w--'); end
            end
            file2 = fullfile(out_group, sprintf('Prop_%s_%s_%s_%s_Fmax%g.%s', Gname, cond, side, pair, FMAX, fmt));
            exportgraphics(fig2, file2, 'Resolution', dpi);

            close([fig1, fig2]);
            fprintf('     ✓ %s | %s : %s (2 figs)\n', cond, side, pair);
        end
    end

    fprintf('  => Images enregistrées dans : %s\n', out_group);
end

fprintf('\n>> Terminé. Dossier racine : %s\n', out_dir);

%% --- utilitaire inline (ternaire string) ---
function s = tern(cond, a, b)
    if cond, s = a; else, s = b; end
end

function key = canon_pair(a,b)
    % Normalise le nom de paire (indépendant de l'ordre)
    u = sort({a,b});               % tri alphabétique
    ab = sprintf('%s_%s', u{1}, u{2});
    % (Optionnel) mapping pour garder tes noms habituels
    switch ab
        case {'GMED_VL','VL_GMED'},             key = 'GMED_VL';
        case {'RF_ST','ST_RF'},                 key = 'RF_ST';
        case {'VL_RF','RF_VL'},                 key = 'VL_RF';
        case {'GM_SOL','SOL_GM'},               key = 'GM_SOL';
        case {'TAprox_TAdist','TAdist_TAprox'}, key = 'TAprox_TAdist';
        otherwise, key = ab; % fallback
    end
end