%% AGRÉGATION DES COHÉRENCES SIGNIFICATIVES + ÉVÉNEMENTS DE MARCHE PAR GROUPE
%    ET VISUALISATION PAR HEATMAPS (VERSION AMÉLIORÉE DU SCRIPT 13)
%
% OBJECTIF :
%   Ce script est une version enrichie du script 13. Il accomplit les mêmes
%   opérations d'agrégation point-à-point des matrices de cohérence
%   significative par groupe d'âge, mais avec deux ajouts majeurs :
%
%   AJOUT 1 — Calcul des événements de marche moyens du groupe :
%     Pendant la boucle d'agrégation, les événements de marche (Toe-Off
%     principal, Toe-Off et Heel Strike controlatéraux) sont extraits depuis
%     les champs GaitEvents_* de chaque fichier Coherence_<PID>.mat. Leur
%     moyenne de groupe est calculée et sauvegardée dans 'GroupEventsAvg'
%     directement dans le fichier .mat de sortie.
%
%   AJOUT 2 — Heatmaps enrichies avec repères d'événements :
%     Les événements moyens du groupe (stockés dans GroupEventsAvg) sont
%     superposés aux heatmaps sous forme de lignes verticales `xline`,
%     avec légende et ticks 0–100% sur l'axe X. Les heatmaps utilisent
%     désormais la colormap `spring` pour la proportion (au lieu de `hot`).
%
% PIPELINE (Partie 1 — Agrégation) :
%   Pour chaque groupe d'âge et chaque fichier Coherence_<PID>.mat :
%     1. Extraction des événements GaitEvents_* (gauche et droite)
%        et accumulation dans EventsRaw.(cond)
%     2. Agrégation point-à-point des matrices CoherenceSignif_<m1>_<m2>
%        (somme, compteur de contributeurs, comptage de présence > 0)
%   Puis, après la boucle de fichiers :
%     3. Construction du côté virtuel 'both' = gauche + droite
%     4. Calcul des moyennes de groupe depuis EventsRaw → GroupEventsAvg
%     5. Sauvegarde : Significant_Coherence_SUM_<Groupe>.mat contenant
%        Significant_Coherence, Significant_Coherence_N, Presence_Count
%        ET GroupEventsAvg (nouveauté vs script 13)
%
% PIPELINE (Partie 2 — Heatmaps) :
%   Pour chaque fichier Significant_Coherence_SUM_<Groupe>.mat :
%     - Heatmap 1 : Somme d'excès de cohérence (colormap parula)
%     - Heatmap 2 : Proportion de sujets avec cohérence > 0 (colormap spring)
%     Les deux heatmaps affichent les lignes d'événements GroupEventsAvg,
%     une légende, et des ticks 0–100% sur l'axe X.
%
% ENTRÉES :
%   - Coherence_<PID>.mat  : Fichiers issus des scripts 7 + 12, contenant :
%       * CoherenceSignif_<m1>_<m2>  : Matrices de cohérence significative
%       * GaitEvents_<m1>_<m2>       : Struct avec main_toeoff,
%                                      opposite_toeoff, opposite_heelstrike
%   - ParticipantGroup.m   : Définit la struct Group (listes de PIDs)
%
% SORTIES :
%   Partie 1 :
%   - Significant_Coherence_SUM_<Groupe>.mat : Un fichier par groupe avec :
%       * Significant_Coherence.(cond).(side).(pair)   : Somme des matrices
%       * Significant_Coherence_N.(cond).(side).(pair) : Nb de contributeurs
%       * Presence_Count.(cond).(side).(pair)           : Nb de sujets > 0
%       * GroupEventsAvg.(cond).main_toeoff             : % moyen du groupe
%       * GroupEventsAvg.(cond).opposite_toeoff         : % moyen du groupe
%       * GroupEventsAvg.(cond).opposite_heelstrike     : % moyen du groupe
%       * Significant_Coherence.meta.Freq               : Axe fréquentiel (Hz)
%     Côtés disponibles : 'left', 'right', 'both' (L+R)
%   Partie 2 :
%   - Sum_<Groupe>_<cond>_both_<pair>_Fmax<n>.png  : Heatmap somme d'excès
%   - Prop_<Groupe>_<cond>_both_<pair>_Fmax<n>.png : Heatmap proportion
%     Avec lignes d'événements, légende et graduation 0–100%
%     Générées dans : Results/<Groupe>/
%
% FONCTIONS LOCALES :
%   - pick_gait_events_block(Str) : Recherche et retourne le premier champ
%                                   GaitEvents_* d'une structure de côté
%   - safe_field_val(ge, field)   : Lecture sécurisée d'un champ scalaire
%                                   (retourne NaN si absent ou non fini)
%
% DIFFÉRENCES AVEC LE SCRIPT 13 :
%   - Extraction et sauvegarde de GroupEventsAvg dans le .mat de sortie
%   - Lignes d'événements sur les heatmaps (xline + légende)
%   - Ticks 0–100% sur l'axe X des heatmaps
%   - Colormap `spring` pour les heatmaps de proportion (vs `hot`)
%   - Presence_Count basé sur M > 0 directement (sans gestion des NaN
%     préalable — différence mineure avec script 13)
% -------------------------------------------------------------------------

%% ===================== PARTIE 1 : AGRÉGATION PAR GROUPE ================
clear; clc; close all;

%% --- Paramètres ---
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
participant_group_file = 'C:\Users\defsil00\Documents\Script\ParticipantGroup.m';
assert(isfolder(all_dir), 'Dossier introuvable: %s', all_dir);

% Groupes à traiter (laisser vide pour tous les groupes définis dans ParticipantGroup.m)
groups_to_run = {}; % ex: {'Adultes','Adolescents'} ou {} pour tous

%% --- Chargement des définitions de groupes ---
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

%% --- Détection des fichiers disponibles ---
files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(files)
    fprintf('Aucun fichier Coherence_*.mat trouvé dans: %s\n', all_dir); return;
end
fprintf('>> %d fichier(s) détecté(s) dans %s\n', numel(files), all_dir);

% Fonction anonyme : extrait le PID depuis le nom de fichier
% ex: 'Coherence_CTL_32.mat' → 'CTL_32'
get_id = @(nm) regexp(nm, '^Coherence_([^\.]+)\.mat$', 'tokens', 'once');

% Prédicat : vrai si le champ est une matrice CoherenceSignif de cycle complet
% Format attendu : 'CoherenceSignif_<m1>_<m2>' (exactement 3 parties)
isFullCycleSignif = @(fn) (startsWith(fn,'CoherenceSignif_') && numel(strsplit(fn,'_'))==3);

sides = {'left','right'};

%% ===================== BOUCLE SUR LES GROUPES ===========================

for ig = 1:numel(groups_to_run)
    Gname      = groups_to_run{ig};
    wanted_ids = upper(string(Group.(Gname))); % IDs du groupe (en majuscules)

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
    Significant_Coherence   = struct(); % Somme point-à-point des matrices
    Significant_Coherence_N = struct(); % Nombre de participants contributeurs (scalaire par paire)
    Presence_Count          = struct(); % Nb de participants avec valeur > 0 en chaque (f,t)
    notes = {};

    % Accumulateur brut des événements de marche, par condition et par participant
    % EventsRaw.(cond).main_toeoff         = [val_p1, val_p2, ...]
    % EventsRaw.(cond).opposite_toeoff     = [...]
    % EventsRaw.(cond).opposite_heelstrike = [...]
    EventsRaw = struct();

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

            % Filtrage des conditions valides (structs avec au moins un côté)
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

                % --- BLOC 1 : Extraction des événements de marche du participant ---
                % Pour chaque côté disponible, récupère les 3 événements clés
                % via pick_gait_events_block (cherche GaitEvents_* dans DATA)
                % puis safe_field_val (lecture sécurisée, NaN si absent/non fini)

                % Côté gauche
                mainL = NaN; oppoToL = NaN; oppoHsL = NaN;
                if isfield(DATA.(condName),'left')
                    geL = pick_gait_events_block(DATA.(condName).left);
                    if ~isempty(geL)
                        mainL   = safe_field_val(geL,'main_toeoff');
                        oppoToL = safe_field_val(geL,'opposite_toeoff');
                        oppoHsL = safe_field_val(geL,'opposite_heelstrike');
                    end
                end

                % Côté droit
                mainR = NaN; oppoToR = NaN; oppoHsR = NaN;
                if isfield(DATA.(condName),'right')
                    geR = pick_gait_events_block(DATA.(condName).right);
                    if ~isempty(geR)
                        mainR   = safe_field_val(geR,'main_toeoff');
                        oppoToR = safe_field_val(geR,'opposite_toeoff');
                        oppoHsR = safe_field_val(geR,'opposite_heelstrike');
                    end
                end

                % Accumulation de la moyenne bilatérale dans EventsRaw
                % (une valeur par participant, unilatéral accepté)
                if ~isfield(EventsRaw, condName), EventsRaw.(condName) = struct(); end

                if any(isfinite([mainL mainR]))
                    if ~isfield(EventsRaw.(condName),'main_toeoff')
                        EventsRaw.(condName).main_toeoff = [];
                    end
                    EventsRaw.(condName).main_toeoff(end+1) = mean([mainL mainR],'omitnan'); %#ok<AGROW>
                end
                if any(isfinite([oppoToL oppoToR]))
                    if ~isfield(EventsRaw.(condName),'opposite_toeoff')
                        EventsRaw.(condName).opposite_toeoff = [];
                    end
                    EventsRaw.(condName).opposite_toeoff(end+1) = mean([oppoToL oppoToR],'omitnan'); %#ok<AGROW>
                end
                if any(isfinite([oppoHsL oppoHsR]))
                    if ~isfield(EventsRaw.(condName),'opposite_heelstrike')
                        EventsRaw.(condName).opposite_heelstrike = [];
                    end
                    EventsRaw.(condName).opposite_heelstrike(end+1) = mean([oppoHsL oppoHsR],'omitnan'); %#ok<AGROW>
                end

                % --- BLOC 2 : Agrégation des matrices de cohérence significative ---
                for s = 1:numel(sides)
                    sideStr = sides{s};
                    if ~isfield(DATA.(condName), sideStr), continue; end
                    Str = DATA.(condName).(sideStr);

                    fns = fieldnames(Str);
                    for iF = 1:numel(fns)
                        fn = fns{iF};
                        if ~isFullCycleSignif(fn), continue; end

                        % Extraction des noms de muscles depuis le nom de champ
                        % ex: 'CoherenceSignif_TAprox_TAdist' → m1='TAprox', m2='TAdist'
                        parts    = strsplit(fn,'_'); % {'CoherenceSignif','m1','m2'}
                        pairName = sprintf('%s_%s', parts{2}, parts{3});

                        M = Str.(fn);
                        if isempty(M) || ndims(M) ~= 2
                            notes{end+1} = sprintf('%s | %s | %s | %s : matrice vide/non-2D', ...
                                files(iFile).name, condName, sideStr, fn);
                            continue;
                        end

                        % Initialisation des nœuds de la structure si nécessaire
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

                        % --- Initialisation (1er participant) ou accumulation ---
                        if ~isfield(Significant_Coherence.(condName).(sideStr), pairName)
                            % Premier participant pour cette paire/condition/côté
                            Significant_Coherence.(condName).(sideStr).(pairName)   = M;
                            Significant_Coherence_N.(condName).(sideStr).(pairName) = 1;
                            % Présence : 1 si la valeur est strictement positive
                            Presence_Count.(condName).(sideStr).(pairName)          = double(M > 0);

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
                            % Participants suivants : vérification dimensionnelle avant addition
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
                                Presence_Count.(condName).(sideStr).(pairName) + double(M > 0);
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

    % --- Construction du côté virtuel 'both' (gauche + droite) ---
    % Pour chaque paire, additionne les matrices gauche et droite ainsi que
    % leurs compteurs respectifs (Significant_Coherence_N et Presence_Count)
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
            pair      = allPairs{ip};
            hasL_pair = hasL && isfield(Significant_Coherence.(cond).left,  pair);
            hasR_pair = hasR && isfield(Significant_Coherence.(cond).right, pair);
            if ~hasL_pair && ~hasR_pair, continue; end

            if hasL_pair && hasR_pair
                % Les deux côtés : vérification dimensionnelle avant addition
                A = Significant_Coherence.(cond).left.(pair);
                B = Significant_Coherence.(cond).right.(pair);
                if ~isequal(size(A), size(B))
                    notes{end+1} = sprintf('%s | both | %s : tailles L/R diff (%s vs %s) -> SKIP', ...
                        cond, pair, mat2str(size(A)), mat2str(size(B)));
                    continue;
                end
                Msum = A + B;

                NL = 0; NR = 0;
                if isfield(Significant_Coherence_N.(cond),'left')  && isfield(Significant_Coherence_N.(cond).left,  pair), NL = Significant_Coherence_N.(cond).left.(pair);  end
                if isfield(Significant_Coherence_N.(cond),'right') && isfield(Significant_Coherence_N.(cond).right, pair), NR = Significant_Coherence_N.(cond).right.(pair); end
                Nsum = NL + NR;

                % Presence_Count both = somme gauche + droite
                if isfield(Presence_Count.(cond),'left')  && isfield(Presence_Count.(cond).left, pair) && ...
                   isfield(Presence_Count.(cond),'right') && isfield(Presence_Count.(cond).right, pair)
                    Psum = Presence_Count.(cond).left.(pair) + Presence_Count.(cond).right.(pair);
                else
                    Psum = [];
                end

            elseif hasL_pair
                % Côté gauche uniquement
                Msum = Significant_Coherence.(cond).left.(pair);
                Nsum = 0;
                if isfield(Significant_Coherence_N.(cond),'left') && isfield(Significant_Coherence_N.(cond).left, pair)
                    Nsum = Significant_Coherence_N.(cond).left.(pair);
                end
                Psum = [];
                if isfield(Presence_Count.(cond),'left') && isfield(Presence_Count.(cond).left, pair)
                    Psum = Presence_Count.(cond).left.(pair);
                end

            else
                % Côté droit uniquement
                Msum = Significant_Coherence.(cond).right.(pair);
                Nsum = 0;
                if isfield(Significant_Coherence_N.(cond),'right') && isfield(Significant_Coherence_N.(cond).right, pair)
                    Nsum = Significant_Coherence_N.(cond).right.(pair);
                end
                Psum = [];
                if isfield(Presence_Count.(cond),'right') && isfield(Presence_Count.(cond).right, pair)
                    Psum = Presence_Count.(cond).right.(pair);
                end
            end

            Significant_Coherence.(cond).both.(pair)   = Msum;
            Significant_Coherence_N.(cond).both.(pair) = Nsum;
            if ~isempty(Psum), Presence_Count.(cond).both.(pair) = Psum; end
        end
    end

    % --- Calcul des événements de marche moyens du groupe (GroupEventsAvg) ---
    % Transforme les vecteurs accumulés dans EventsRaw en scalaires moyens.
    % GroupEventsAvg est ensuite sauvegardé dans le .mat de sortie pour être
    % utilisé directement par la Partie 2 sans relecture des fichiers individuels.
    GroupEventsAvg = struct();
    evConds = fieldnames(EventsRaw);
    for ic = 1:numel(evConds)
        cond = evConds{ic};
        GroupEventsAvg.(cond) = struct();

        % Toe-Off principal moyen du groupe
        if isfield(EventsRaw.(cond),'main_toeoff') && ~isempty(EventsRaw.(cond).main_toeoff)
            GroupEventsAvg.(cond).main_toeoff = mean(EventsRaw.(cond).main_toeoff, 'omitnan');
        else
            GroupEventsAvg.(cond).main_toeoff = NaN;
        end

        % Toe-Off controlatéral moyen du groupe
        if isfield(EventsRaw.(cond),'opposite_toeoff') && ~isempty(EventsRaw.(cond).opposite_toeoff)
            GroupEventsAvg.(cond).opposite_toeoff = mean(EventsRaw.(cond).opposite_toeoff, 'omitnan');
        else
            GroupEventsAvg.(cond).opposite_toeoff = NaN;
        end

        % Heel Strike controlatéral moyen du groupe
        if isfield(EventsRaw.(cond),'opposite_heelstrike') && ~isempty(EventsRaw.(cond).opposite_heelstrike)
            GroupEventsAvg.(cond).opposite_heelstrike = mean(EventsRaw.(cond).opposite_heelstrike, 'omitnan');
        else
            GroupEventsAvg.(cond).opposite_heelstrike = NaN;
        end

        fprintf('   [%s] main_TO=%.2f%% | oppo_TO=%.2f%% | oppo_HS=%.2f%%\n', cond, ...
            GroupEventsAvg.(cond).main_toeoff, ...
            GroupEventsAvg.(cond).opposite_toeoff, ...
            GroupEventsAvg.(cond).opposite_heelstrike);
    end

    % --- Sauvegarde du fichier de synthèse (avec GroupEventsAvg) ---
    % Note : GroupEventsAvg est une variable supplémentaire par rapport au
    % script 13, ce qui rend ce fichier .mat directement exploitable
    % pour les visualisations sans relire les données individuelles.
    out_file = fullfile(all_dir, sprintf('Significant_Coherence_SUM_%s.mat', Gname));
    save(out_file, 'Significant_Coherence', 'Significant_Coherence_N', ...
                   'Presence_Count', 'GroupEventsAvg', '-v7.3');
    fprintf('>> Sauvegardé: %s\n', out_file);
end

fprintf('\n== Terminé pour groupes: %s ==\n', strjoin(groups_to_run, ', '));

%% ===================== PARTIE 2 : HEATMAPS PAR GROUPE ==================
% Génère deux heatmaps par combinaison groupe/condition/paire, enrichies
% par rapport au script 13 avec :
%   - Lignes verticales des événements de marche (xline)
%   - Légende des événements
%   - Ticks 0–100% sur l'axe X
clear; clc; close all;

% --- Paramètres ---
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
out_dir = 'C:\Users\defsil00\Documents\Script\Results';
side    = 'both';   % Côté à visualiser (gauche + droite combinés)
FMAX    = 80;       % Fréquence maximale affichée dans les heatmaps (Hz)
dpi     = 300;      % Résolution d'export en PNG
fmt     = 'png';    % Format : 'png' | 'tiff' | 'jpg'

% Restriction optionnelle à certains groupes (vide = tous les fichiers trouvés)
groups_to_plot = {};

% --- Détection des fichiers de synthèse par groupe ---
gfiles = dir(fullfile(all_dir, 'Significant_Coherence_SUM_*.mat'));
% Exclut l'éventuel fichier global sans suffixe de groupe
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
    if ~isfield(S,'Significant_Coherence') || ...
       ~isfield(S,'Significant_Coherence_N') || ...
       ~isfield(S,'Presence_Count')
        fprintf('  [skip] Structures manquantes dans %s\n', gfiles(ig).name); continue;
    end
    Significant_Coherence   = S.Significant_Coherence;
    Significant_Coherence_N = S.Significant_Coherence_N;
    Presence_Count          = S.Presence_Count;

    % Chargement de GroupEventsAvg (avec compatibilité fichiers anciens sans ce champ)
    if isfield(S,'GroupEventsAvg')
        GroupEventsAvg = S.GroupEventsAvg;
    else
        GroupEventsAvg = struct(); % Fichier produit par script 13 : pas d'événements
    end

    % Axe fréquentiel (Hz)
    hasHz = isfield(Significant_Coherence,'meta') && ...
            isfield(Significant_Coherence.meta,'Freq');
    if hasHz
        f_all = Significant_Coherence.meta.Freq(:);
    else
        warning('(%s) meta.Freq manquant : axe Y en indices de bins.', Gname);
    end

    % Création du dossier de sortie propre au groupe
    out_group = fullfile(out_dir, Gname);
    if ~exist(out_group,'dir'), mkdir(out_group); end

    % --- Boucle conditions → paires ---
    condList = setdiff(fieldnames(Significant_Coherence), {'meta'});
    for ic = 1:numel(condList)
        cond = condList{ic};
        if ~isfield(Significant_Coherence.(cond), side), continue; end

        pairList = fieldnames(Significant_Coherence.(cond).(side));
        fprintf('  %s : %d paires\n', cond, numel(pairList));

        for ip = 1:numel(pairList)
            pair = pairList{ip};

            % Vérification de la disponibilité de toutes les structures
            if ~isfield(Significant_Coherence_N.(cond), side) || ...
               ~isfield(Significant_Coherence_N.(cond).(side), pair) || ...
               ~isfield(Presence_Count.(cond), side) || ...
               ~isfield(Presence_Count.(cond).(side), pair)
                fprintf('   [skip] Données incomplètes: %s | %s | %s\n', cond, side, pair);
                continue;
            end

            SumMap = Significant_Coherence.(cond).(side).(pair);   % Somme des excès [nFreq × 1000]
            N      = Significant_Coherence_N.(cond).(side).(pair); % Nb de contributeurs (scalaire)
            P      = Presence_Count.(cond).(side).(pair);           % Nb de sujets > 0 par (f,t)

            [nFreq, nTime] = size(SumMap);
            t = linspace(0, 100, nTime); % Axe temporel normalisé (% du cycle)

            % Construction de l'axe fréquentiel et tri croissant
            if hasHz && numel(f_all) == nFreq
                f = f_all(:);
            else
                f = (1:nFreq)'; % Fallback : index de bin
            end
            [fc, idx] = sort(f, 'ascend');
            SumMap_c  = SumMap(idx, :);
            P_c       = P(idx, :);

            % Découpage à FMAX
            mask = hasHz && (fc <= FMAX);
            if ~any(mask), mask = true(size(fc)); end
            fvis = fc(mask);
            SumV = SumMap_c(mask, :);
            Pv   = P_c(mask, :);
            Prop = Pv ./ max(N, 1); % Proportion [0..1] : Nb sujets actifs / Nb total

            % --- Fonction locale pour les lignes d'événements ---
            % Factorisée ici en commentaires pour éviter la duplication de code :
            % Pour chaque heatmap, on trace 3 xline si les valeurs sont finies :
            %   main_toeoff       : --k (noir tirets)      | Toe-Off principal
            %   opposite_toeoff   : :   (magenta pointillé) | TO controlatéral
            %   opposite_heelstrike: -. (cyan dash-point)   | HS controlatéral

            % --- Heatmap 1 : Somme d'excès de cohérence ---
            fig1 = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
            imagesc(t, fvis, SumV);
            set(gca,'YDir','normal'); axis tight; colorbar;
            xlabel('Gait cycle (%)');
            if hasHz, ylabel('Frequency (Hz)'); else, ylabel('Frequency (bin)'); end
            title(sprintf('%s | %s | %s — Sum of excess — %s', ...
                cond, side, strrep(pair,'_','-'), FMAX, Gname));
            colormap(parula);
            hold on;

            % Lignes horizontales délimitant les bandes fréquentielles (Alpha/Beta/Gamma)
            if hasHz
                yb = [8 12 13 30 31 60];
                yb = yb(yb >= min(fvis) & yb <= max(fvis));
                if ~isempty(yb), yline(yb,'w--'); end
            end

            % Graduation de l'axe X (0 à 100% par pas de 10%)
            xticks(0:10:100);
            xticklabels(arrayfun(@(x)sprintf('%d%%',x), 0:10:100, 'UniformOutput',false));

            % Superposition des lignes d'événements de marche du groupe
            legItems = []; legLabels = {};
            if isfield(GroupEventsAvg, cond)
                E = GroupEventsAvg.(cond);
                % TO principal (trait noir tirets)
                if isfield(E,'main_toeoff') && isfinite(E.main_toeoff)
                    h = xline(E.main_toeoff,'--k','LineWidth',1.5,'DisplayName','main toe-off');
                    legItems(end+1) = h; legLabels{end+1} = 'main toe-off';
                end
                % TO controlatéral (pointillé magenta)
                if isfield(E,'opposite_toeoff') && isfinite(E.opposite_toeoff)
                    h = xline(E.opposite_toeoff,':','Color',[0.85 0 0.85],'LineWidth',1.5,'DisplayName','opposite toe-off');
                    legItems(end+1) = h; legLabels{end+1} = 'opposite toe-off';
                end
                % HS controlatéral (dash-point cyan)
                if isfield(E,'opposite_heelstrike') && isfinite(E.opposite_heelstrike)
                    h = xline(E.opposite_heelstrike,'-.','Color',[0 0.6 1],'LineWidth',1.5,'DisplayName','opposite heel-strike');
                    legItems(end+1) = h; legLabels{end+1} = 'opposite heel-strike';
                end
            end
            if ~isempty(legItems)
                legend(legItems, legLabels, 'Location','northoutside','Orientation','horizontal');
            end

            file1 = fullfile(out_group, sprintf('Sum_%s_%s_%s_%s_Fmax%g.%s', ...
                Gname, cond, side, pair, FMAX, fmt));
            exportgraphics(fig1, file1, 'Resolution', dpi);

            % --- Heatmap 2 : Proportion de sujets avec cohérence > 0 ---
            % Colormap 'spring' (vs 'hot' dans script 13) — gamme rose/magenta
            fig2 = figure('Visible','off','Color','w','Position',[100 100 1200 500]);
            imagesc(t, fvis, Prop);
            set(gca,'YDir','normal'); axis tight; colorbar; caxis([0 1]);
            xlabel('Gait cycle (%)');
            if hasHz, ylabel('Frequency (Hz)'); else, ylabel('Frequency (bin)'); end
            title(sprintf('%s | %s | %s — Proportion — %s', ...
                cond, side, strrep(pair,'_','-'), FMAX, Gname));
            colormap(spring);
            hold on;

            % Mêmes lignes de bandes fréquentielles
            if hasHz
                yb = [8 12 13 30 31 60];
                yb = yb(yb >= min(fvis) & yb <= max(fvis));
                if ~isempty(yb), yline(yb,'w--'); end
            end

            % Même graduation 0–100%
            xticks(0:10:100);
            xticklabels(arrayfun(@(x)sprintf('%d%%',x), 0:10:100, 'UniformOutput',false));

            % Mêmes lignes d'événements (identiques à la heatmap 1)
            legItems = []; legLabels = {};
            if isfield(GroupEventsAvg, cond)
                E = GroupEventsAvg.(cond);
                if isfield(E,'main_toeoff') && isfinite(E.main_toeoff)
                    h = xline(E.main_toeoff,'--k','LineWidth',1.5,'DisplayName','main toe-off');
                    legItems(end+1) = h; legLabels{end+1} = 'main toe-off';
                end
                if isfield(E,'opposite_toeoff') && isfinite(E.opposite_toeoff)
                    h = xline(E.opposite_toeoff,':','Color',[0.85 0 0.85],'LineWidth',1.5,'DisplayName','opposite toe-off');
                    legItems(end+1) = h; legLabels{end+1} = 'opposite toe-off';
                end
                if isfield(E,'opposite_heelstrike') && isfinite(E.opposite_heelstrike)
                    h = xline(E.opposite_heelstrike,'-.','Color',[0 0.6 1],'LineWidth',1.5,'DisplayName','opposite heel-strike');
                    legItems(end+1) = h; legLabels{end+1} = 'opposite heel-strike';
                end
            end
            if ~isempty(legItems)
                legend(legItems, legLabels, 'Location','northoutside','Orientation','horizontal');
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

function ge = pick_gait_events_block(Str)
% PICK_GAIT_EVENTS_BLOCK  Retourne la structure GaitEvents_* d'un côté.
%
%   Cherche le premier champ commençant par 'GaitEvents_' dans la structure
%   du côté (DATA.(cond).(side)). En l'absence d'un tel champ, tente le
%   champ de référence 'GaitEvents_TAprox_TAdist' (paire la plus souvent
%   présente). Retourne [] si rien n'est trouvé ou si le contenu n'est pas
%   une struct.
%
%   Entrée : Str — DATA.(cond).(side), structure d'un côté
%   Sortie : ge  — struct GaitEvents, ou [] si introuvable

    ge = [];
    if ~isstruct(Str), return; end
    fns = fieldnames(Str);
    idx = find(startsWith(fns,'GaitEvents_'), 1, 'first');
    if ~isempty(idx)
        ge = Str.(fns{idx});
    elseif isfield(Str,'GaitEvents_TAprox_TAdist')
        ge = Str.('GaitEvents_TAprox_TAdist'); 
    end
    if ~isstruct(ge), ge = []; end
end

function v = safe_field_val(ge, fieldname)
% SAFE_FIELD_VAL  Lecture sécurisée d'un champ scalaire fini dans une struct.
%
%   Retourne la valeur de ge.(fieldname) uniquement si :
%     - ge est une struct
%     - le champ fieldname existe
%     - la valeur est un scalaire numérique fini (non NaN, non Inf)
%   Dans tous les autres cas, retourne NaN.
%
%   Entrées : ge        — struct GaitEvents (ou tout autre struct)
%             fieldname — nom du champ à lire (string)
%   Sortie  : v         — scalaire double, ou NaN

    v = NaN;
    if isstruct(ge) && isfield(ge, fieldname)
        val = ge.(fieldname);
        if isscalar(val) && isfinite(val)
            v = val;
        end
    end
end