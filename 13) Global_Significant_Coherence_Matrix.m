%% AGRÉGATION POINT-À-POINT DES COHÉRENCES SIGNIFICATIVES PAR GROUPE D'ÂGE
%
% OBJECTIF :
%   Ce script agrège, pour chaque groupe d'âge, les matrices de cohérence
%   significative (CoherenceSignif_<m1>_<m2>, produites par le script 12)
%   de tous les participants du groupe. Pour chaque combinaison
%   condition/côté/paire, il calcule :
%     - La somme point-à-point des matrices de cohérence significative
%     - Le nombre de contributeurs (participants ayant fourni cette matrice)
%     - Un comptage de présence : nb de participants avec cohérence > 0 en
%       chaque point (f,t) — utilisé pour calculer la proportion de sujets
%
%   Un côté virtuel 'both' est ensuite construit par addition des côtés
%   gauche et droit.
%
%   La deuxième partie du script charge ces matrices agrégées et génère des
%   heatmaps de la somme d'excès et de la proportion de sujets, pour chaque
%   combinaison groupe/condition/paire, exportées en PNG haute résolution.
%
% ENTRÉES :
%   - Coherence_<PID>.mat   : Fichiers de résultats du script 12, contenant
%                             les champs CoherenceSignif_<m1>_<m2> dans DATA.
%   - ParticipantGroup.m    : Script définissant la struct Group avec les
%                             listes de participants par groupe d'âge.
%
% SORTIES :
%   Partie 1 :
%   - Significant_Coherence_SUM_<Groupe>.mat : Un fichier par groupe, contenant :
%       * Significant_Coherence.(cond).(side).(pair)   : Somme des matrices
%       * Significant_Coherence_N.(cond).(side).(pair) : Nb de contributeurs
%       * Presence_Count.(cond).(side).(pair)           : Nb de sujets > 0
%                                                         en chaque point (f,t)
%       * Significant_Coherence.meta.Freq              : Axe fréquentiel (Hz)
%     Les trois côtés sont disponibles : 'left', 'right', 'both' (L+R)
%   Partie 2 :
%   - Sum_<Groupe>_<cond>_<side>_<pair>_Fmax<n>.<fmt>  : Heatmap de la somme
%   - Prop_<Groupe>_<cond>_<side>_<pair>_Fmax<n>.<fmt> : Heatmap de proportion
%     Générées dans : Visualisation_Coherence/<Groupe>/
%
% PARAMÈTRES CONFIGURABLES :
%   - groups_to_run  : Liste de groupes à traiter (vide = tous)
%   - groups_to_plot : Liste de groupes à visualiser (vide = tous les .mat trouvés)
%   - FMAX           : Fréquence maximale affichée dans les heatmaps (Hz)
%   - make_backup / dry_run : Options de sécurité (non utilisées ici)
%
% FONCTIONS LOCALES :
%   - canon_pair(a,b) : Normalise l'ordre des muscles dans la clé de paire
%   - tern(cond,a,b)  : Opérateur ternaire pour les labels d'axes
% -------------------------------------------------------------------------

%% ===================== PARTIE 1 : AGRÉGATION PAR GROUPE ================
clear; clc; close all;

%% --- Paramètres ---
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
participant_group_file = 'C:\Users\defsil00\Documents\Script\ParticipantGroup.m';
assert(isfolder(all_dir), 'Dossier introuvable: %s', all_dir);

% Groupes à traiter (laisser vide pour traiter tous les groupes définis)
groups_to_run = {}; % ex: {'JeunesEnfants'} ou {} pour tous

% --- Chargement des définitions de groupes ---
assert(exist(participant_group_file,'file')==2, ...
    'Fichier ParticipantGroup.m introuvable: %s', participant_group_file);
run(participant_group_file); % → crée la variable struct Group
assert(exist('Group','var')==1 && isstruct(Group), ...
    'ParticipantGroup.m doit définir une variable struct "Group".');

all_groups = fieldnames(Group);
if isempty(groups_to_run)
    groups_to_run = all_groups;
else
    % Suppression des groupes demandés mais absents de la définition
    miss = setdiff(groups_to_run, all_groups);
    if ~isempty(miss)
        warning('Groupes inconnus ignorés: %s', strjoin(miss, ', '));
        groups_to_run = setdiff(groups_to_run, miss);
    end
end
if isempty(groups_to_run)
    warning('Aucun groupe valide à traiter.'); return;
end

% --- Détection des fichiers disponibles ---
files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(files)
    fprintf('Aucun fichier Coherence_*.mat trouvé dans: %s\n', all_dir); return;
end
fprintf('>> %d fichier(s) détecté(s) dans %s\n', numel(files), all_dir);

% Fonctions anonymes utilitaires
% Extrait le PID depuis le nom de fichier (ex: 'Coherence_CTL_04.mat' → 'CTL_04')
get_id = @(nm) regexp(nm, '^Coherence_([^\.]+)\.mat$', 'tokens', 'once');

% Prédicat : vrai si un champ est une matrice CoherenceSignif de cycle complet
% Format attendu : 'CoherenceSignif_<m1>_<m2>' (exactement 3 parties)
isFullCycleSignif = @(fn) (startsWith(fn,'CoherenceSignif_') && numel(strsplit(fn,'_'))==3);

sides = {'left','right'};

%% ===================== BOUCLE SUR LES GROUPES ===========================

for ig = 1:numel(groups_to_run)
    Gname      = groups_to_run{ig};
    wanted_ids = upper(string(Group.(Gname))); % IDs du groupe (en majuscules pour la comparaison)

    % --- Sélection des fichiers appartenant au groupe ---
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
        % Signalement des participants attendus mais absents du dossier
        missing_ids = setdiff(wanted_ids, unique(file_ids));
        if ~isempty(missing_ids)
            fprintf('   (absents dans le dossier: %s)\n', strjoin(cellstr(missing_ids), ', '));
        end
    end

    % --- Initialisation des structures d'agrégation pour ce groupe ---
    % Significant_Coherence   : somme point-à-point des matrices
    % Significant_Coherence_N : nombre de participants contributeurs (scalaire par paire)
    % Presence_Count          : nombre de participants avec valeur > 0 en chaque (f,t)
    Significant_Coherence   = struct();
    Significant_Coherence_N = struct();
    Presence_Count          = struct();
    notes = {};

    % ===================== BOUCLE SUR LES FICHIERS DU GROUPE ============
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

            % Filtrage des conditions valides
            condNames = fieldnames(DATA);
            mask = false(size(condNames));
            for k = 1:numel(condNames)
                nm = condNames{k};
                mask(k) = isstruct(DATA.(nm)) && ...
                          (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
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

                        % Extraction des noms de muscles et normalisation de la clé
                        % canon_pair garantit un ordre alphabétique constant
                        % (indépendant de l'ordre dans le champ)
                        parts    = strsplit(fn,'_'); % {'CoherenceSignif','m1','m2'}
                        pairName = canon_pair(parts{2}, parts{3});

                        M = Str.(fn);
                        if isempty(M) || ndims(M) ~= 2
                            notes{end+1} = sprintf('%s | %s | %s | %s : matrice vide/non-2D', ...
                                files(iFile).name, condName, sideStr, fn);
                            continue;
                        end

                        % --- Préparation de la matrice avant agrégation ---
                        if ~isfloat(M), M = double(M); end

                        % Masque de présence : points avec cohérence > 0 et non-NaN
                        % (utilisé pour le comptage de sujets actifs par point)
                        mask_pos = (~isnan(M)) & (M > 0);

                        % Remplacement des NaN par 0 pour ne pas polluer la somme
                        M(isnan(M)) = 0;

                        % --- Initialisation des nœuds si nécessaire ---
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

                        % --- Initialisation ou accumulation ---
                        if ~isfield(Significant_Coherence.(condName).(sideStr), pairName)
                            % Premier participant pour cette paire/condition/côté
                            Significant_Coherence.(condName).(sideStr).(pairName)   = M;
                            Significant_Coherence_N.(condName).(sideStr).(pairName) = 1;
                            Presence_Count.(condName).(sideStr).(pairName)          = double(mask_pos);

                            % Stockage de l'axe fréquentiel (une seule fois, si disponible)
                            if ~isfield(Significant_Coherence,'meta') || ...
                               ~isfield(Significant_Coherence.meta,'Freq')
                                if isfield(Str,'Freq') && ~isempty(Str.Freq) && ...
                                   numel(Str.Freq) == size(M,1)
                                    Significant_Coherence.meta.Freq = Str.Freq(:);
                                end
                            end
                            fprintf('     [+] init %s | %s | %s (%s)\n', ...
                                condName, sideStr, pairName, mat2str(size(M)));
                        else
                            % Participants suivants : vérification de la compatibilité
                            % dimensionnelle avant addition
                            Msum = Significant_Coherence.(condName).(sideStr).(pairName);
                            if ~isequal(size(Msum), size(M))
                                notes{end+1} = sprintf( ...
                                    '%s | %s | %s | %s : taille différente (%s vs %s) -> SKIP', ...
                                    files(iFile).name, condName, sideStr, pairName, ...
                                    mat2str(size(Msum)), mat2str(size(M)));
                                continue;
                            end
                            % Accumulation de la somme, du compteur et de la présence
                            Significant_Coherence.(condName).(sideStr).(pairName) = Msum + M;
                            Significant_Coherence_N.(condName).(sideStr).(pairName) = ...
                                Significant_Coherence_N.(condName).(sideStr).(pairName) + 1;
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
    end % Fichiers du groupe

    % --- Récapitulatif par groupe ---
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
        if numel(notes) > max_show
            fprintf('  ... (+%d supplémentaires)\n', numel(notes)-max_show);
        end
    end

    % --- Construction du côté virtuel 'both' (gauche + droit) ---
    % Pour chaque paire, additionne les matrices gauche et droite
    % ainsi que leurs compteurs de participants et de présence
    condList = setdiff(fieldnames(Significant_Coherence), {'meta'});
    for iC = 1:numel(condList)
        cond = condList{iC};
        hasL = isfield(Significant_Coherence.(cond),'left');
        hasR = isfield(Significant_Coherence.(cond),'right');
        if ~hasL && ~hasR, continue; end

        % Union des paires disponibles sur l'un ou l'autre côté
        Lpairs = {}; Rpairs = {};
        if hasL, Lpairs = fieldnames(Significant_Coherence.(cond).left);  end
        if hasR, Rpairs = fieldnames(Significant_Coherence.(cond).right); end
        allPairs = unique([Lpairs; Rpairs]);

        % Initialisation des nœuds 'both'
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
            pair     = allPairs{ip};
            hasL_pair = hasL && isfield(Significant_Coherence.(cond).left,  pair);
            hasR_pair = hasR && isfield(Significant_Coherence.(cond).right, pair);
            if ~hasL_pair && ~hasR_pair, continue; end

            if hasL_pair && hasR_pair
                % Les deux côtés disponibles : vérification dimensionnelle avant addition
                A = Significant_Coherence.(cond).left.(pair);
                B = Significant_Coherence.(cond).right.(pair);
                if ~isequal(size(A), size(B))
                    notes{end+1} = sprintf('%s | both | %s : tailles L/R diff (%s vs %s) -> SKIP', ...
                        cond, pair, mat2str(size(A)), mat2str(size(B)));
                    continue;
                end
                A(isnan(A)) = 0; B(isnan(B)) = 0; % Sécurisation finale
                Msum = A + B;

                NL = 0; NR = 0;
                if isfield(Significant_Coherence_N.(cond),'left')  && isfield(Significant_Coherence_N.(cond).left,  pair), NL = Significant_Coherence_N.(cond).left.(pair);  end
                if isfield(Significant_Coherence_N.(cond),'right') && isfield(Significant_Coherence_N.(cond).right, pair), NR = Significant_Coherence_N.(cond).right.(pair); end
                Nsum = NL + NR;

                % Addition des comptages de présence si disponibles sur les deux côtés
                if isfield(Presence_Count.(cond),'left')  && isfield(Presence_Count.(cond).left, pair) && ...
                   isfield(Presence_Count.(cond),'right') && isfield(Presence_Count.(cond).right, pair)
                    Pl = Presence_Count.(cond).left.(pair);  Pl(isnan(Pl)) = 0;
                    Pr = Presence_Count.(cond).right.(pair); Pr(isnan(Pr)) = 0;
                    Psum = Pl + Pr;
                else
                    Psum = [];
                end

            elseif hasL_pair
                % Côté gauche uniquement
                Msum = Significant_Coherence.(cond).left.(pair); Msum(isnan(Msum)) = 0;
                Nsum = 0;
                if isfield(Significant_Coherence_N.(cond),'left') && isfield(Significant_Coherence_N.(cond).left, pair)
                    Nsum = Significant_Coherence_N.(cond).left.(pair);
                end
                Psum = [];
                if isfield(Presence_Count.(cond),'left') && isfield(Presence_Count.(cond).left, pair)
                    Psum = Presence_Count.(cond).left.(pair); Psum(isnan(Psum)) = 0;
                end

            else
                % Côté droit uniquement
                Msum = Significant_Coherence.(cond).right.(pair); Msum(isnan(Msum)) = 0;
                Nsum = 0;
                if isfield(Significant_Coherence_N.(cond),'right') && isfield(Significant_Coherence_N.(cond).right, pair)
                    Nsum = Significant_Coherence_N.(cond).right.(pair);
                end
                Psum = [];
                if isfield(Presence_Count.(cond),'right') && isfield(Presence_Count.(cond).right, pair)
                    Psum = Presence_Count.(cond).right.(pair); Psum(isnan(Psum)) = 0;
                end
            end

            % Enregistrement du côté 'both'
            Significant_Coherence.(cond).both.(pair)   = Msum;
            Significant_Coherence_N.(cond).both.(pair) = Nsum;
            if ~isempty(Psum), Presence_Count.(cond).both.(pair) = Psum; end
        end
    end

    % --- Sauvegarde du fichier de synthèse pour ce groupe ---
    out_file = fullfile(all_dir, sprintf('Significant_Coherence_SUM_%s.mat', Gname));
    save(out_file, 'Significant_Coherence', 'Significant_Coherence_N', 'Presence_Count', '-v7.3');
    fprintf('>> Sauvegardé: %s\n', out_file);
end

fprintf('\n== Terminé pour groupes: %s ==\n', strjoin(groups_to_run, ', '));

%% ===================== PARTIE 2 : HEATMAPS PAR GROUPE ==================
% Génère deux heatmaps par combinaison groupe/condition/paire :
%   1. Somme d'excès de cohérence (SumMap) — colormap parula
%   2. Proportion de sujets avec cohérence > 0 (Prop = P / N) — colormap hot
clear; clc; close all;

% --- Paramètres ---
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
out_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\Visualisation_Coherence';
side    = 'both'; % Côté à visualiser ('left', 'right', ou 'both')
FMAX    = 80;     % Fréquence maximale affichée dans les heatmaps (Hz)
dpi     = 300;
fmt     = 'png';  % Format d'export : 'png' | 'tiff' | 'jpg'

% Restriction optionnelle à certains groupes (vide = tous les fichiers trouvés)
groups_to_plot = {};

% --- Détection des fichiers de synthèse par groupe ---
gfiles = dir(fullfile(all_dir, 'Significant_Coherence_SUM_*.mat'));
% Exclut l'éventuel fichier non-groupé (sans suffixe de groupe)
gfiles = gfiles(~strcmp({gfiles.name}, 'Significant_Coherence_SUM.mat'));

if ~isempty(groups_to_plot)
    keep = false(1, numel(gfiles));
    for i = 1:numel(gfiles)
        tok = regexp(gfiles(i).name, '^Significant_Coherence_SUM_(.+)\.mat$', 'tokens','once');
        if ~isempty(tok)
            keep(i) = any(strcmpi(tok{1}, groups_to_plot));
        end
    end
    gfiles = gfiles(keep);
end
if isempty(gfiles)
    error('Aucun fichier groupe trouvé (pattern Significant_Coherence_SUM_*.mat) dans %s', all_dir);
end

if ~exist(out_dir,'dir'), mkdir(out_dir); end

% ===================== BOUCLE DE VISUALISATION ==========================
for ig = 1:numel(gfiles)
    tok = regexp(gfiles(ig).name, '^Significant_Coherence_SUM_(.+)\.mat$', 'tokens','once');
    if isempty(tok)
        fprintf('[skip] Fichier inattendu: %s\n', gfiles(ig).name); continue;
    end
    Gname = tok{1};
    fprintf('\n=== GROUPE: %s ===\n', Gname);

    % Chargement des structures d'agrégation
    S = load(fullfile(gfiles(ig).folder, gfiles(ig).name));
    if ~isfield(S, 'Significant_Coherence') || ...
       ~isfield(S, 'Significant_Coherence_N') || ...
       ~isfield(S, 'Presence_Count')
        fprintf('  [skip] Structures manquantes dans %s\n', gfiles(ig).name); continue;
    end
    Significant_Coherence   = S.Significant_Coherence;
    Significant_Coherence_N = S.Significant_Coherence_N;
    Presence_Count          = S.Presence_Count;

    % Axe fréquentiel (Hz) — utilisé pour l'axe Y et le découpage à FMAX
    hasHz = isfield(Significant_Coherence,'meta') && ...
            isfield(Significant_Coherence.meta,'Freq');
    if hasHz
        f_all = Significant_Coherence.meta.Freq(:);
    else
        warning('(%s) meta.Freq manquant : axe Y en indices de bins.', Gname);
    end

    % Création du dossier de sortie propre au groupe
    out_group = fullfile(out_dir, Gname);
    if ~exist(out_group, 'dir'), mkdir(out_group); end

    % --- Boucle conditions → paires ---
    condList = setdiff(fieldnames(Significant_Coherence), {'meta'});
    for ic = 1:numel(condList)
        cond = condList{ic};
        if ~isfield(Significant_Coherence.(cond), side), continue; end

        pairList = fieldnames(Significant_Coherence.(cond).(side));
        fprintf('  %s : %d paires\n', cond, numel(pairList));

        for ip = 1:numel(pairList)
            pair = pairList{ip};

            % Vérification de la disponibilité de toutes les structures nécessaires
            if ~isfield(Significant_Coherence_N.(cond), side) || ...
               ~isfield(Significant_Coherence_N.(cond).(side), pair) || ...
               ~isfield(Presence_Count.(cond), side) || ...
               ~isfield(Presence_Count.(cond).(side), pair)
                fprintf('   [skip] Données incomplètes: %s | %s | %s\n', cond, side, pair);
                continue;
            end

            SumMap = Significant_Coherence.(cond).(side).(pair); % Somme des excès [nFreq × 1000]
            SumMap(isnan(SumMap)) = 0; % Protection finale
            N      = Significant_Coherence_N.(cond).(side).(pair); % Nb de contributeurs (scalaire)
            P      = Presence_Count.(cond).(side).(pair);           % Nb de sujets > 0 par (f,t)

            [nFreq, nTime] = size(SumMap);
            t = linspace(0, 100, nTime); % Axe temporel normalisé (% du cycle)

            % Construction de l'axe fréquentiel
            if hasHz && numel(f_all) == nFreq
                f = f_all(:);
            else
                f = (1:nFreq)'; % Fallback : index de bin
            end

            % Tri par fréquence croissante et découpage à FMAX
            [fc, idx] = sort(f, 'ascend');
            SumMap_c  = SumMap(idx, :);
            P_c       = P(idx, :);
            mask      = hasHz && (fc <= FMAX); % Si pas de Hz : tout afficher
            if ~any(mask), mask = true(size(fc)); end

            fvis = fc(mask);
            SumV = SumMap_c(mask, :);
            Pv   = P_c(mask, :);
            Prop = Pv ./ max(N, 1); % Proportion [0..1] : Nb sujets actifs / Nb total

            % --- Heatmap 1 : Somme d'excès de cohérence ---
            fig1 = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
            SumV(isnan(SumV)) = 0; % Sécurité ultime
            imagesc(t, fvis, SumV);
            set(gca,'YDir','normal'); axis tight; colorbar;
            xlabel('Gait cycle (%)');
            ylabel(tern(hasHz, 'Frequency (Hz)', 'Frequency (bin)'));
            title(sprintf('%s | %s | %s — Sum of excess — %s', ...
                cond, side, strrep(pair,'_','-'), Gname));
            colormap(parula);
            % Tracé des lignes de délimitation des bandes fréquentielles
            if hasHz
                yb = [8 12 13 30 31 60];
                yb = yb(yb >= min(fvis) & yb <= max(fvis));
                if ~isempty(yb), hold on; yline(yb,'w--'); end
            end
            file1 = fullfile(out_group, sprintf('Sum_%s_%s_%s_%s_Fmax%g.%s', ...
                Gname, cond, side, pair, FMAX, fmt));
            exportgraphics(fig1, file1, 'Resolution', dpi);

            % --- Heatmap 2 : Proportion de sujets avec cohérence > 0 ---
            fig2 = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
            imagesc(t, fvis, Prop);
            set(gca,'YDir','normal'); axis tight; colorbar; caxis([0 1]);
            xlabel('Gait cycle (%)');
            ylabel(tern(hasHz, 'Frequency (Hz)', 'Frequency (bin)'));
            title(sprintf('%s | %s | %s — Proportion — %s', ...
                cond, side, strrep(pair,'_','-'), Gname));
            colormap(hot);
            if hasHz
                yb = [8 12 13 30 31 60];
                yb = yb(yb >= min(fvis) & yb <= max(fvis));
                if ~isempty(yb), hold on; yline(yb,'w--'); end
            end
            file2 = fullfile(out_group, sprintf('Prop_%s_%s_%s_%s_Fmax%g.%s', ...
                Gname, cond, side, pair, FMAX, fmt));
            exportgraphics(fig2, file2, 'Resolution', dpi);

            close([fig1, fig2]);
            fprintf('     ✓ %s | %s : %s (2 figs)\n', cond, side, pair);
        end
    end
    fprintf('  => Images enregistrées dans : %s\n', out_group);
end

fprintf('\n>> Terminé. Dossier racine : %s\n', out_dir);

%% ===================== FONCTIONS LOCALES ================================

function s = tern(cond, a, b)
% TERN  Opérateur ternaire : retourne a si cond est vrai, b sinon.
%   Usage : s = tern(isValid, 'label_vrai', 'label_faux')
    if cond, s = a; else, s = b; end
end

function key = canon_pair(a, b)
% CANON_PAIR  Normalise la clé d'une paire de muscles (indépendant de l'ordre).
%
%   Trie les deux noms alphabétiquement puis applique un mapping optionnel
%   pour conserver les noms habituels du projet (ex: 'VL_RF' plutôt que 'RF_VL').
%
%   Entrées : a, b — noms des deux muscles (strings)
%   Sortie  : key  — clé normalisée de la forme 'muscle1_muscle2'

    u  = sort({a, b}); % Tri alphabétique pour un ordre canonique
    ab = sprintf('%s_%s', u{1}, u{2});

    % Mapping vers les noms de paires conventionnels du projet
    switch ab
        case {'GMED_VL','VL_GMED'},             key = 'GMED_VL';
        case {'RF_ST','ST_RF'},                  key = 'RF_ST';
        case {'VL_RF','RF_VL'},                  key = 'VL_RF';
        case {'GM_SOL','SOL_GM'},                key = 'GM_SOL';
        case {'TAprox_TAdist','TAdist_TAprox'},  key = 'TAprox_TAdist';
        otherwise, key = ab; % Fallback : ordre alphabétique
    end
end