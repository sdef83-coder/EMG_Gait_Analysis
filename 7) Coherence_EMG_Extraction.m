%% ANALYSE DE COHÉRENCE EMG-EMG PAR TRANSFORMÉE EN ONDELETTES CONTINUES
%
% OBJECTIF :
%   Ce script calcule la cohérence temps-fréquence entre paires de muscles
%   ipsilatéraux (même jambe) à partir de signaux EMG bruts issus de fichiers
%   .c3d. L'analyse repose sur la transformée en ondelettes continues (Morlet)
%   et produit des matrices de cohérence découpées selon les sous-phases du
%   cycle de marche (Loading Response, Mid-Stance, Pre-Swing, Swing).
%
%   Pour garantir la comparabilité inter-conditions, le nombre de cycles
%   utilisés est égalisé au minimum disponible parmi les trois conditions
%   (Plat, Medium, High), après exclusion des cycles artefactés (outliers).
%   La sélection des cycles est randomisée de façon reproductible (graine RNG
%   fixée) pour assurer la réplicabilité des résultats.
%
% ENTRÉES :
%   - Fichiers .c3d      : Signaux EMG et événements de marche (Heel Strike /
%                          Toe Off) par essai, condition et participant.
%   - *_MATRIX.mat       : Fichier pré-calculé par participant contenant la
%                          matrice de signaux et la structure CYCLES_OUTLIERS
%                          (cycles à exclure).
%   - Association.m      : Script de mapping liant les labels musculaires
%                          (ex: 'EMG_TAprox') aux capteurs analogiques réels
%                          du système Vicon, séparément pour la jambe gauche
%                          (sensor_association_left) et droite
%                          (sensor_association_right).
%   - Mapping-EMG.xlsm   : Fichier Excel contenant la validité (code = 2)
%                          de chaque canal EMG par participant et condition,
%                          utilisé pour filtrer les paires invalides.
%   - Fonctions utilitaires :
%       indiceLeft / indiceRight     : Détection des Heel Strikes (début de cycle)
%       indiceLeftTO / indiceRightTO : Détection des Toe-Offs
%       WaveletParameters            : Paramétrage de la transformée ondelette
%       wavelet                      : Calcul de la TFR (Transformée en
%                                      Fréquences et Résolutions)
%       filtrage                     : Filtre passe-bande appliqué aux EMG
%
% SORTIES :
%   - Coherence_<PID>.mat : Fichier .mat par participant contenant la
%                           structure DATA avec, pour chaque paire de muscles,
%                           condition et côté :
%       * Matrices de cohérence complète et par sous-phase
%       * Moyennes de cohérence par bande fréquentielle (Alpha/Beta/Gamma)
%       * Spectres de puissance individuels et inter-spectre
%       * Indices des sous-phases (Loading, MidStance, PreSwing, Swing)
%       * Événements de marche moyens (TO/HS) utilisés pour le découpage
%       * Seuil de significativité statistique de la cohérence
%       * Métadonnées (nb de cycles, graine RNG, cible d'égalisation)
%   - Figures .png        : Une figure par paire/condition/côté montrant les
%                           spectres de puissance, la cohérence significative
%                           et les profils fréquentiels par sous-phase.
%
% DÉPENDANCES :
%   - BTK Toolbox (btkReadAcquisition, btkGetAnalogs, btkGetAnalogFrequency)
%   - Fonctions maison dans le dossier 'coherence_analysis'
% -------------------------------------------------------------------------

clear; close all; clc;

%% ===================== CHEMINS & DÉPENDANCES ============================

% Ajout des fonctions d'analyse de cohérence au path MATLAB
addpath(genpath('C:\Users\defsil00\Documents\Script\Functions\coherence_analysis'));

% Répertoire de travail contenant les fichiers .c3d du groupe d'âge à analyser
cd('C:\Users\defsil00\Documents\Script\Data\jeunes_enfants'); % <- adapter selon la catégorie d'âge

% Lecture du fichier Excel de mapping EMG (validité des canaux par participant)
% 'raw' contient les données brutes (texte + numérique) nécessaires aux vérifications
[num, txt, raw] = xlsread('C:\Users\defsil00\Documents\Script\Mapping-EMG.xlsm');

%% ===================== PARAMÉTRAGE UTILISATEUR ==========================

% Identifiant(s) du ou des participant(s) à traiter
Participant = {['CTL_63']};   % Exemple : un seul participant

% Conditions expérimentales et leur décalage de colonne dans le fichier Excel
% de mapping (colonne de base = colonne du participant dans 'raw')
%   Plat   → décalage +1
%   Medium → décalage +3
%   High   → décalage +5
Condition = {
    'Plat',   1;
    'Medium', 3;
    'High',   5
};

% Numéros d'essais à inclure dans l'analyse
Essai = {'01','02','03','04', '05', '06', '07', '08', '09', '10'}; 

% Indices de ligne dans 'raw' pour chaque muscle
% (permettent de vérifier la validité du canal dans le fichier Excel)
raw_indices = struct( ...
    'TAprox', 13, ...
    'TAdist', 14, ...
    'SOL',    15, ...
    'GM',     16, ...
    'VL',     17, ...
    'RF',     18, ...
    'ST',     19, ...
    'GMED',   20 ...
);

% Seuil de significativité pour le test de cohérence (formule de Rosenberg)
% seuil = 1 - alpha^(1 / (N_cycles - 1))
alpha = 0.05;

%% ===================== PARAMÈTRES DE L'ANALYSE ONDELETTE ================

% Plage fréquentielle analysée (Hz) et résolution spectrale
FreqMin    = 1;
FreqMax    = 400;
Resolution = 1;

% Numéro d'onde de l'ondelette de Morlet (compromis temps-fréquence)
WaveNumber = 7;

% Définition des bandes fréquentielles d'intérêt
FreqBands = struct( ...
    'Alpha', [8  12], ...
    'Beta',  [13  30], ...
    'Gamma', [31  60] ...
);

%% ===================== REPRODUCTIBILITÉ — GRAINE RNG ====================

% Fixation de la graine du générateur aléatoire pour que la sélection
% aléatoire des cycles (égalisation inter-conditions) soit reproductible
rng(0,'twister');
seedState = rng;
fprintf('RNG seed: %d (%s)\n', seedState.Seed, seedState.Type);

%% ===================== BOUCLE PRINCIPALE PAR PARTICIPANT ================

for iP = 1:numel(Participant)

    % --- Chargement des données pré-calculées du participant ---
    % Le fichier *_MATRIX.mat contient la matrice EMG normalisée et la
    % structure CYCLES_OUTLIERS (cycles artefactés à exclure)
    load(['C:\Users\defsil00\Documents\Script\ORIGINALS\' ...
          Participant{iP} '_MATRIX.mat']);

    % Chargement du mapping capteurs ↔ muscles pour ce participant
    % Produit : sensor_association_left et sensor_association_right
    run('C:\Users\defsil00\Documents\Script\Association.m');

    % Initialisation de la structure des outliers si absente du fichier .mat
    if ~exist('CYCLES_OUTLIERS','var')
        CYCLES_OUTLIERS = struct();
    end

    % Création du dossier de sortie pour les résultats du participant
    output_dir = fullfile('C:\Users\defsil00\Documents\Script\Results\Coherence', Participant{iP});
    if ~exist(output_dir, 'dir'); mkdir(output_dir); end

    % --- Définition des paires de muscles à analyser ---
    % Chaque ligne = une paire {muscle1, muscle2} dont on calculera la cohérence
    muscle_pairs = {
        'TAprox', 'TAdist';  % Tibialis Anterior proximal - distal
        'VL',     'RF';      % Vastus Lateralis - Rectus Femoris
        'GM',     'SOL';     % Gastrocnemius Medialis - Soleus
        'GMED',   'RF';      % Gluteus Medius - Rectus Femoris
        'GMED',   'VL';      % Gluteus Medius - Vastus Lateralis
        'RF',     'ST'       % Rectus Femoris - Semitendinosus
        % ... ajouter d'autres paires si nécessaire
    };

    % --- Construction de la table des paires valides (Pairs) ---
    % Pour chaque paire de muscles et chaque côté (gauche/droit), vérifie
    % que les deux capteurs existent dans le mapping du participant.
    % Pairs = {sensor1, sensor2, nom_m1, nom_m2, côté, row_m1, row_m2}
    Pairs = {};
    for i = 1:size(muscle_pairs, 1)
        m1 = muscle_pairs{i,1};
        m2 = muscle_pairs{i,2};

        for side = {'left','right'}
            side_str = side{1};

            % Sélection du mapping capteur selon le côté
            if strcmp(side_str,'left')
                assoc = sensor_association_left;
            else
                assoc = sensor_association_right;
            end

            % Ajout à la table si les deux capteurs et leurs indices sont disponibles
            if isfield(assoc, ['EMG_' m1]) && isfield(assoc, ['EMG_' m2]) && ...
               isfield(raw_indices, m1) && isfield(raw_indices, m2)

                Pairs(end+1,:) = { ...
                    assoc.(['EMG_' m1]), ...   % Nom du canal analogique muscle 1
                    assoc.(['EMG_' m2]), ...   % Nom du canal analogique muscle 2
                    m1, ...                    % Label muscle 1
                    m2, ...                    % Label muscle 2
                    side_str, ...              % Côté ('left' ou 'right')
                    raw_indices.(m1), ...      % Ligne dans 'raw' pour muscle 1
                    raw_indices.(m2) ...       % Ligne dans 'raw' pour muscle 2
                };
            end
        end
    end

    % ===================== CALCUL DES ÉVÉNEMENTS MOYENS (TO/HS) =========
    % Pour chaque condition et chaque côté, calcule les pourcentages moyens
    % du cycle auxquels surviennent les événements clés :
    %   - Toe-Off de la jambe principale       (fin de phase d'appui)
    %   - Toe-Off de la jambe opposée          (début double appui controlatéral)
    %   - Heel Strike de la jambe opposée      (reprise d'appui controlatéral)
    % Ces événements seront utilisés pour délimiter les sous-phases du cycle.
    GaitEvents = struct();

    for iC = 1:size(Condition,1)
        for main_side = {'left','right'}
            main_side_str = main_side{1};

            % Détermination du côté opposé
            if strcmp(main_side_str,'left')
                opposite_side_str = 'right';
            else
                opposite_side_str = 'left';
            end

            % Accumulateurs pour les pourcentages d'événements sur tous les essais
            main_toe_offs_all         = [];
            opposite_toe_offs_all     = [];
            opposite_heel_strikes_all = [];

            for iEs = 1:numel(Essai)
                file = [Participant{iP} '_' Condition{iC,1} '_' Essai{iEs} '.c3d'];
                if ~isfile(file); continue; end

                data    = btkReadAcquisition(file);
                analogs = btkGetAnalogs(data);
                FreqS   = btkGetAnalogFrequency(data);
                FreqVicon = 100; % Fréquence d'acquisition vidéo Vicon (Hz)

                % Signal EMG temporaire (requis comme argument par les fonctions
                % d'indexation, mais non utilisé dans le calcul des événements)
                flds     = fieldnames(analogs);
                temp_emg = analogs.(flds{1});
                EMG_temp = [temp_emg, temp_emg]; % Matrice factice à 2 colonnes

                % Récupération des indices d'événements — jambe principale
                if strcmp(main_side_str,'left')
                    main_heel_strikes = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_temp);
                    main_toe_offs     = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                else
                    main_heel_strikes = indiceRight(data, analogs, FreqS, FreqVicon, EMG_temp);
                    main_toe_offs     = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                end

                % Récupération des indices d'événements — jambe opposée
                if strcmp(opposite_side_str,'left')
                    opposite_heel_strikes = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_temp);
                    opposite_toe_offs     = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                else
                    opposite_heel_strikes = indiceRight(data, analogs, FreqS, FreqVicon, EMG_temp);
                    opposite_toe_offs     = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                end

                % --- Boucle sur les cycles de l'essai ---
                % Pour chaque cycle valide (HS → HS suivant contenant un TO),
                % calcule les pourcentages d'occurrence des événements dans le cycle
                for ii = 1:min(numel(main_heel_strikes)-1, numel(main_toe_offs))

                    % Vérifie que le TO se situe bien entre les deux HS consécutifs
                    if main_toe_offs(ii) > main_heel_strikes(ii) && main_toe_offs(ii) < main_heel_strikes(ii+1)

                        cycle_start  = main_heel_strikes(ii);
                        cycle_end    = main_heel_strikes(ii+1);
                        cycle_length = cycle_end - cycle_start;

                        % TO jambe principale (% du cycle)
                        main_toeoff_percent = ((main_toe_offs(ii) - cycle_start) / cycle_length) * 100;

                        % TO jambe opposée : 1er TO controlatéral dans la fenêtre du cycle
                        opposite_to_in_cycle = opposite_toe_offs( ...
                            opposite_toe_offs > cycle_start & opposite_toe_offs < cycle_end);
                        if ~isempty(opposite_to_in_cycle)
                            opposite_toeoff_percent = ((opposite_to_in_cycle(1) - cycle_start) / cycle_length) * 100;
                        else
                            opposite_toeoff_percent = NaN;
                        end

                        % HS jambe opposée : 1er HS controlatéral dans la fenêtre du cycle
                        opposite_hs_in_cycle = opposite_heel_strikes( ...
                            opposite_heel_strikes > cycle_start & opposite_heel_strikes < cycle_end);
                        if ~isempty(opposite_hs_in_cycle)
                            opposite_heelstrike_percent = ((opposite_hs_in_cycle(1) - cycle_start) / cycle_length) * 100;
                        else
                            opposite_heelstrike_percent = NaN;
                        end

                        % Accumulation uniquement si les valeurs sont dans la plage [0, 100]
                        if main_toeoff_percent > 0 && main_toeoff_percent < 100
                            main_toe_offs_all(end+1) = main_toeoff_percent; 
                        end
                        if ~isnan(opposite_toeoff_percent) && opposite_toeoff_percent > 0 && opposite_toeoff_percent < 100
                            opposite_toe_offs_all(end+1) = opposite_toeoff_percent; 
                        end
                        if ~isnan(opposite_heelstrike_percent) && opposite_heelstrike_percent > 0 && opposite_heelstrike_percent < 100
                            opposite_heel_strikes_all(end+1) = opposite_heelstrike_percent;
                        end
                    end
                end
            end % Essais

            % Calcul des moyennes (ou valeurs par défaut si aucun événement trouvé)
            % Valeurs de référence bibliographiques utilisées par défaut :
            %   TO principal ≈ 60%, TO controlatéral ≈ 10%, HS controlatéral ≈ 50%
            if ~isempty(main_toe_offs_all)
                GaitEvents.(Condition{iC,1}).(main_side_str).main_toeoff = mean(main_toe_offs_all);
            else
                GaitEvents.(Condition{iC,1}).(main_side_str).main_toeoff = 60;
            end

            if ~isempty(opposite_toe_offs_all)
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_toeoff = mean(opposite_toe_offs_all);
            else
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_toeoff = 10;
            end

            if ~isempty(opposite_heel_strikes_all)
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_heelstrike = mean(opposite_heel_strikes_all);
            else
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_heelstrike = 50;
            end

            % Affichage des événements moyens calculés
            fprintf('Événements moyens calculés pour %s - %s:\n', Condition{iC,1}, main_side_str);
            fprintf('  Main TO : %.1f%% | Opp TO : %.1f%% | Opp HS : %.1f%%\n', ...
                GaitEvents.(Condition{iC,1}).(main_side_str).main_toeoff, ...
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_toeoff, ...
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_heelstrike);
        end
    end

    % ===================== TRAITEMENT PRINCIPAL : COHÉRENCE =============
    % Boucle sur chaque paire de muscles.
    % Pour chaque paire, le traitement suit 4 étapes :
    %   1) Pré-calcul des indices de cycles valides par condition et côté
    %   2) Détermination de la cible d'égalisation (min des totaux valides)
    %   3) Construction des masques de sélection aléatoire (KeepMask)
    %   4) Calcul des spectres et de la cohérence sur les cycles retenus
    DATA = struct();

    for iPair = 1:size(Pairs, 1)
        m1 = Pairs{iPair,3};
        m2 = Pairs{iPair,4};
        pair_label = [m1 '_' m2];

        % Recherche de la colonne du participant dans le tableau Excel 'raw'
        ii_col = [];
        for jj = 1:size(raw,2)
            if strcmp(raw{1,jj}, Participant{iP})
                ii_col = jj;
                break;
            end
        end
        if isempty(ii_col)
            warning('Participant %s introuvable dans "raw". Paire %s ignorée.', Participant{iP}, pair_label);
            continue;
        end

        % --- ÉTAPE 1 : Pré-calcul des cycles valides par condition et côté ---
        % Pour chaque condition et côté, identifie les indices de cycles dont
        % aucun des deux muscles de la paire n'est marqué comme outlier.
        ValidIdx     = struct(); % Indices des cycles valides (vecteur par cond/côté)
        NtotalSide   = struct(); % Nombre total de cycles (valides + invalides) par cond/côté
        L_total_cond = zeros(size(Condition,1),1); % Nb total de cycles valides par condition (gauche + droit)

        for iC = 1:size(Condition,1)
            condName = Condition{iC,1};

            for side_str_cell = {'left','right'}
                s = side_str_cell{1};
                side_offset = strcmp(s,'right'); % 0 = gauche, 1 = droit (décalage de colonne)

                % Vérifie que les deux muscles sont valides (code = 2 dans le mapping Excel)
                is_ok = false;
                if strcmp(s,'left')
                    assoc_side = sensor_association_left;
                else
                    assoc_side = sensor_association_right;
                end

                if isfield(assoc_side, ['EMG_' m1]) && isfield(assoc_side, ['EMG_' m2]) && ...
                   isfield(raw_indices, m1) && isfield(raw_indices, m2)
                    try
                        if raw{Pairs{iPair,6}, ii_col + Condition{iC,2} + side_offset} == 2 && ...
                           raw{Pairs{iPair,7}, ii_col + Condition{iC,2} + side_offset} == 2
                            is_ok = true;
                        end
                    catch
                        is_ok = false;
                    end
                end

                % Si la paire est valide, identifie les cycles sans outlier
                % (tutu_vec = 0 signifie que les deux muscles sont propres sur ce cycle)
                if is_ok && isfield(CYCLES_OUTLIERS, Participant{iP}) && ...
                   isfield(CYCLES_OUTLIERS.(Participant{iP}), condName) && ...
                   isfield(CYCLES_OUTLIERS.(Participant{iP}).(condName), s)

                    BadCycles = CYCLES_OUTLIERS.(Participant{iP}).(condName).(s);
                    tutu_vec  = BadCycles.(['EMG_' m1]) + BadCycles.(['EMG_' m2]); % 0 = bon sur les 2 muscles
                    ValidIdx.(condName).(s) = find(tutu_vec==0);
                    NtotalSide.(condName).(s) = numel(tutu_vec);
                else
                    ValidIdx.(condName).(s) = [];
                    NtotalSide.(condName).(s) = 0;
                end
            end

            % Nombre total de cycles valides pour cette condition (gauche + droit combinés)
            L_total_cond(iC) = numel(ValidIdx.(condName).left) + numel(ValidIdx.(condName).right);
        end

        % --- ÉTAPE 2 : Détermination de la cible d'égalisation ---
        % La cible = minimum du nombre de cycles valides (gauche+droit) parmi
        % toutes les conditions strictement non-nulles. Cette approche garantit
        % que le même nombre de cycles est utilisé pour toutes les conditions,
        % évitant ainsi un biais lié aux différences de taille d'échantillon.
        pos_mask = L_total_cond > 0;
        if any(pos_mask)
            L_target_eff = min(L_total_cond(pos_mask)); % Minimum strictement positif
        else
            L_target_eff = 0; % Aucune condition n'a de cycles valides
        end

        % Avertissement uniquement si AUCUNE condition ne dispose de cycles
        if L_target_eff == 0
            fprintf('⚠ Paire %s : aucune condition n''a de cycles valides (%s). Paire ignorée.\n', pair_label, Participant{iP});
        end

        % Sauvegarde de la cible dans la structure DATA pour traçabilité
        for iiC = 1:size(Condition,1)
            condName_ii = Condition{iiC,1};
            DATA.(condName_ii).meta.(['L_target_' m1 '_' m2]) = L_target_eff;
            if isfield(NtotalSide, condName_ii)
                DATA.(condName_ii).meta.(['Ltot_' m1 '_' m2 '_left'])  = NtotalSide.(condName_ii).left;
                DATA.(condName_ii).meta.(['Ltot_' m1 '_' m2 '_right']) = NtotalSide.(condName_ii).right;
            end
        end

        % --- ÉTAPE 3 : Construction des masques de sélection (KeepMask) ---
        % KeepMask indique, pour chaque cycle et chaque condition, s'il doit
        % être inclus dans l'analyse. Si le nombre de cycles valides dépasse
        % la cible, une sélection aléatoire est effectuée (randomisée via
        % la graine RNG fixée ci-dessus pour la reproductibilité).
        KeepMask = struct();
        for iC = 1:size(Condition,1)
            condName = Condition{iC,1};
            nL = NtotalSide.(condName).left;
            nR = NtotalSide.(condName).right;

            % Initialisation : aucun cycle retenu par défaut
            KeepMask.(condName).left  = false(nL,1);
            KeepMask.(condName).right = false(nR,1);

            idxL = ValidIdx.(condName).left;
            idxR = ValidIdx.(condName).right;

            % Table unifiée gauche + droit pour la sélection aléatoire conjointe
            union_sides = [
                [repmat({'left'},  numel(idxL), 1), num2cell(idxL(:))];
                [repmat({'right'}, numel(idxR), 1), num2cell(idxR(:))]
            ];
            n_union = size(union_sides,1);

            if n_union == 0
                continue; % Aucun cycle valide pour cette condition
            end

            if L_target_eff == 0
                % Cas dégénéré : on conserve tous les cycles valides disponibles
                if ~isempty(idxL), KeepMask.(condName).left(idxL)  = true; end
                if ~isempty(idxR), KeepMask.(condName).right(idxR) = true; end
            else
                if n_union <= L_target_eff
                    % Pas besoin de sous-échantillonner : on garde tout
                    if ~isempty(idxL), KeepMask.(condName).left(idxL)  = true; end
                    if ~isempty(idxR), KeepMask.(condName).right(idxR) = true; end
                else
                    % Tirage aléatoire sans remise pour atteindre L_target_eff cycles
                    sel    = randperm(n_union, L_target_eff);
                    chosen = union_sides(sel,:);
                    for k = 1:size(chosen,1)
                        if strcmp(chosen{k,1}, 'left')
                            KeepMask.(condName).left(chosen{k,2}) = true;
                        else
                            KeepMask.(condName).right(chosen{k,2}) = true;
                        end
                    end
                end
            end
        end

        % --- ÉTAPE 4 : Calcul des spectres et de la cohérence ---
        % Pour chaque condition, lit les fichiers .c3d, applique la TFR ondelette
        % sur les cycles sélectionnés par KeepMask, accumule les spectres de
        % puissance et l'inter-spectre, puis calcule la cohérence normalisée.
        for iC = 1:size(Condition,1)
            condName = Condition{iC,1};

            % Côté de la paire (défini lors de la construction de la table Pairs)
            side_str = Pairs{iPair,5};
            Side = {side_str, double(strcmp(side_str,'right'))}; % {'left',0} ou {'right',1}

            % Vérifie à nouveau la validité du canal dans le mapping Excel
            is_ok = false;
            try
                if raw{Pairs{iPair,6}, ii_col + Condition{iC,2} + Side{1,2}} == 2 && ...
                   raw{Pairs{iPair,7}, ii_col + Condition{iC,2} + Side{1,2}} == 2
                    is_ok = true;
                end
            catch
                is_ok = false;
            end
            if ~is_ok, continue; end

            % Récupération du vecteur de marquage des outliers pour CE côté
            % (tutu_vec(k) = 0 → cycle k est propre sur les deux muscles)
            if isfield(CYCLES_OUTLIERS, Participant{iP}) && ...
               isfield(CYCLES_OUTLIERS.(Participant{iP}), condName) && ...
               isfield(CYCLES_OUTLIERS.(Participant{iP}).(condName), side_str)
                BadCycles = CYCLES_OUTLIERS.(Participant{iP}).(condName).(side_str);
                tutu_vec  = BadCycles.(['EMG_' m1]) + BadCycles.(['EMG_' m2]);
            else
                tutu_vec  = [];
            end

            % Initialisation des accumulateurs spectraux
            % (dimensionnés dynamiquement lors du 1er cycle retenu)
            PowSpec_s1     = [];
            PowSpec_s2     = [];
            cross_spectrum = [];
            NcycleP = 0; % Compteur global de cycles (tous essais confondus)
            NcycleC = 0; % Compteur de cycles effectivement conservés

            for iEs = 1:numel(Essai)
                file = [Participant{iP} '_' condName '_' Essai{iEs} '.c3d'];
                if ~isfile(file), continue; end

                data    = btkReadAcquisition(file);
                analogs = btkGetAnalogs(data);

                % Lecture des deux signaux EMG de la paire
                EMG = [];
                EMG(:,1) = analogs.(Pairs{iPair,1}); % Muscle 1
                EMG(:,2) = analogs.(Pairs{iPair,2}); % Muscle 2

                FreqS     = btkGetAnalogFrequency(data);
                FreqVicon = 100;

                % Détection des indices de début de cycle (Heel Strikes)
                if strcmp(Side{1},'left')
                    Cycles = indiceLeft(data, analogs, FreqS, FreqVicon, EMG);
                else
                    Cycles = indiceRight(data, analogs, FreqS, FreqVicon, EMG);
                end

                % Paramétrage de la transformée en ondelettes continues
                Args = WaveletParameters(FreqMin, FreqMax, Resolution, WaveNumber, FreqS);

                % Prétraitements du signal EMG
                EMG = filtrage(EMG, FreqS, 8, 400); % Filtre passe-bande [8–400 Hz]
                EMG = EMG - mean(EMG,1);             % Centrage (suppression de la composante continue)

                % Calcul de la Transformée en Fréquences et Résolutions (TFR)
                % pour chaque muscle (ondelette de Morlet)
                [TFR(:,:,1), period, ~, ~] = wavelet(EMG(:,1), Args.DT, Args.Pad, Args.DJ, Args.S0, Args.J1, Args.Mother, Args.Cycles);
                [TFR(:,:,2), ~,     ~, ~]  = wavelet(EMG(:,2), Args.DT, Args.Pad, Args.DJ, Args.S0, Args.J1, Args.Mother, Args.Cycles);
                Freq = 1 ./ period; % Conversion période → fréquence (Hz)

                % --- Boucle sur les cycles de l'essai ---
                for iCycles = 1:(numel(Cycles)-1)
                    NcycleP = NcycleP + 1; % Incrémentation de l'index global

                    % Double condition d'inclusion :
                    %   1) Le cycle est valide (non marqué outlier sur les 2 muscles)
                    %   2) Le cycle est sélectionné par le masque d'égalisation
                    is_valid = (NcycleP > numel(tutu_vec)) || (tutu_vec(NcycleP)==0);
                    in_mask  = (NcycleP <= numel(KeepMask.(condName).(side_str))) && ...
                                KeepMask.(condName).(side_str)(NcycleP);

                    if is_valid && in_mask

                        % Découpe de la TFR sur la fenêtre du cycle
                        TFR_cycle(:,:,1) = TFR(:, Cycles(iCycles):Cycles(iCycles+1), 1);
                        TFR_cycle(:,:,2) = TFR(:, Cycles(iCycles):Cycles(iCycles+1), 2);

                        % Interpolation à 1000 points (normalisation de la longueur du cycle)
                        % → permet de comparer des cycles de durées différentes
                        [X, Y]   = meshgrid(1:size(TFR_cycle,2), 1:size(TFR_cycle,1));
                        [Xq, Yq] = meshgrid(linspace(1,size(TFR_cycle,2),1000), 1:size(TFR_cycle,1));

                        TFR_int(:,:,1) = interp2(X, Y, TFR_cycle(:,:,1), Xq, Yq, 'spline');
                        TFR_int(:,:,2) = interp2(X, Y, TFR_cycle(:,:,2), Xq, Yq, 'spline');

                        % Initialisation des accumulateurs au premier cycle retenu
                        if isempty(PowSpec_s1)
                            PowSpec_s1     = zeros(size(TFR_int,1), size(TFR_int,2));
                            PowSpec_s2     = zeros(size(TFR_int,1), size(TFR_int,2));
                            cross_spectrum = zeros(size(TFR_int,1), size(TFR_int,2));
                        end

                        % Accumulation des sommes (division par N faite après la boucle)
                        PowSpec_s1     = PowSpec_s1     + abs(TFR_int(:,:,1)).^2;          % Σ|TFR1|²
                        PowSpec_s2     = PowSpec_s2     + abs(TFR_int(:,:,2)).^2;          % Σ|TFR2|²
                        cross_spectrum = cross_spectrum + (TFR_int(:,:,1)) .* conj(TFR_int(:,:,2)); % ΣTFr1·TFR2*

                        NcycleC = NcycleC + 1;

                        clear TFR_int TFR_cycle
                    end
                end % Cycles

                clear TFR
            end % Essais

            % Si aucun cycle n'a été retenu, passe à la paire/condition suivante
            if NcycleC == 0
                warning('Aucun cycle retenu après égalisation | %s-%s | %s | %s', m1, m2, condName, side_str);
                continue;
            end

            % --- Calcul de la cohérence ---
            % Cohérence = |Σ(TFR1·TFR2*)|² / (Σ|TFR1|² · Σ|TFR2|²)
            % (estimateur de cohérence par sommation de Welch adapté aux ondelettes)
            Pxx_sum = PowSpec_s1;
            Pyy_sum = PowSpec_s2;
            Pxy_sum = cross_spectrum;
            L_cycles = NcycleC;

            % Moyennage sur le nombre de cycles retenus
            PowSpec_s1_mean     = PowSpec_s1 / NcycleC;
            PowSpec_s2_mean     = PowSpec_s2 / NcycleC;
            cross_spectrum_mean = cross_spectrum / NcycleC;

            % Matrice de cohérence normalisée (valeurs entre 0 et 1)
            Coherence = abs(cross_spectrum_mean).^2 ./ (PowSpec_s1_mean .* PowSpec_s2_mean);

            % --- Délimitation des sous-phases du cycle ---
            % Conversion des pourcentages d'événements en indices sur 1000 points
            ge = GaitEvents.(condName).(side_str);

            main_toeoff_index         = max(1,   min(round(ge.main_toeoff      * 10), 1000));
            opposite_toeoff_index     = max(1,   min(round(ge.opposite_toeoff  * 10), max(1, main_toeoff_index-1)));
            opposite_heelstrike_index = max(opposite_toeoff_index+1, ...
                                            min(round(ge.opposite_heelstrike * 10), max(2, main_toeoff_index-1)));

            % Définition des sous-phases (indices sur le vecteur normalisé 1000 pts)
            %   Loading Response : HS principal → TO controlatéral
            %   Mid-Stance       : TO controlatéral → HS controlatéral
            %   Pre-Swing        : HS controlatéral → TO principal
            %   Swing            : TO principal → HS principal suivant
            loading_response_indices = 1:opposite_toeoff_index;
            midstance_indices        = (opposite_toeoff_index+1):opposite_heelstrike_index;
            preswing_indices         = (opposite_heelstrike_index+1):main_toeoff_index;
            swing_phase_indices      = (main_toeoff_index+1):1000;

            % --- Moyennes de cohérence par bande fréquentielle et sous-phase ---
            mean_coherence_loading   = struct();
            mean_coherence_midstance = struct();
            mean_coherence_preswing  = struct();
            mean_coherence_swing     = struct();
            mean_coherence_full      = struct();

            for bandName = fieldnames(FreqBands)'
                band      = bandName{1};
                range     = FreqBands.(band);
                idx_band  = find(Freq >= range(1) & Freq <= range(2));

                % Moyenne sur tout le cycle
                mean_coherence_full.(band) = mean(mean(Coherence(idx_band,:), 2), 'omitnan');

                % Moyennes par sous-phase (NaN si la sous-phase est vide)
                if ~isempty(loading_response_indices)
                    mean_coherence_loading.(band) = mean(mean(Coherence(idx_band, loading_response_indices), 2), 'omitnan');
                else
                    mean_coherence_loading.(band) = NaN;
                end

                if ~isempty(midstance_indices)
                    mean_coherence_midstance.(band) = mean(mean(Coherence(idx_band, midstance_indices), 2), 'omitnan');
                else
                    mean_coherence_midstance.(band) = NaN;
                end

                if ~isempty(preswing_indices)
                    mean_coherence_preswing.(band) = mean(mean(Coherence(idx_band, preswing_indices), 2), 'omitnan');
                else
                    mean_coherence_preswing.(band) = NaN;
                end

                if ~isempty(swing_phase_indices)
                    mean_coherence_swing.(band) = mean(mean(Coherence(idx_band, swing_phase_indices), 2), 'omitnan');
                else
                    mean_coherence_swing.(band) = NaN;
                end

                % Enregistrement des scalaires dans DATA
                DATA.(condName).(side_str).(['MeanCoherence_' band '_' m1 '_' m2])                 = mean_coherence_full.(band);
                DATA.(condName).(side_str).(['MeanCoherence_LoadingResponse_' band '_' m1 '_' m2]) = mean_coherence_loading.(band);
                DATA.(condName).(side_str).(['MeanCoherence_MidStance_' band '_' m1 '_' m2])       = mean_coherence_midstance.(band);
                DATA.(condName).(side_str).(['MeanCoherence_PreSwing_' band '_' m1 '_' m2])        = mean_coherence_preswing.(band);
                DATA.(condName).(side_str).(['MeanCoherence_Swing_' band '_' m1 '_' m2])           = mean_coherence_swing.(band);
            end

            % --- Seuil de significativité de la cohérence (Rosenberg et al.) ---
            % Formule : seuil = 1 - alpha^(1 / (N - 1))
            % Indéfini si N ≤ 1 (un seul cycle ne suffit pas à estimer la cohérence)
            if NcycleC <= 1
                seuil = NaN;
            else
                seuil = 1 - alpha^(1/(NcycleC - 1));
            end

            % --- Sauvegarde des résultats dans la structure DATA ---

            % Événements de marche utilisés pour le découpage en sous-phases
            DATA.(condName).(side_str).(['GaitEvents_' m1 '_' m2]) = struct( ...
                'main_toeoff',         ge.main_toeoff, ...
                'opposite_toeoff',     ge.opposite_toeoff, ...
                'opposite_heelstrike', ge.opposite_heelstrike);

            % Matrices de cohérence découpées par sous-phase
            DATA.(condName).(side_str).(['Coherence_LoadingResponse_' m1 '_' m2]) = Coherence(:, loading_response_indices);
            DATA.(condName).(side_str).(['Coherence_MidStance_'       m1 '_' m2]) = Coherence(:, midstance_indices);
            DATA.(condName).(side_str).(['Coherence_PreSwing_'        m1 '_' m2]) = Coherence(:, preswing_indices);
            DATA.(condName).(side_str).(['Coherence_Swing_'           m1 '_' m2]) = Coherence(:, swing_phase_indices);

            % Indices (en points sur 1000) délimitant chaque sous-phase
            DATA.(condName).(side_str).(['LoadingResponseIndices_' m1 '_' m2]) = loading_response_indices;
            DATA.(condName).(side_str).(['MidStanceIndices_'       m1 '_' m2]) = midstance_indices;
            DATA.(condName).(side_str).(['PreSwingIndices_'        m1 '_' m2]) = preswing_indices;
            DATA.(condName).(side_str).(['SwingIndices_'           m1 '_' m2]) = swing_phase_indices;

            % Spectres de puissance individuels et croisé (sommes sur tous les cycles)
            DATA.(condName).(side_str).(['PowSpec' m1]) = PowSpec_s1;
            DATA.(condName).(side_str).(['PowSpec' m2]) = PowSpec_s2;
            DATA.(condName).(side_str).(['Cross_Spectrum_' m1 '_' m2]) = abs(cross_spectrum / NcycleC).^2; % |<TFR1·TFR2*>|²

            % Matrice de cohérence complète (toutes fréquences × 1000 points)
            DATA.(condName).(side_str).(['Coherence_' m1 '_' m2]) = Coherence;

            % Nombre de cycles effectivement utilisés pour ce calcul
            DATA.(condName).(side_str).(['Ncycle_' m1 '_' m2]) = NcycleC;

            % Métadonnées d'égalisation et de reproductibilité
            DATA.(condName).meta.(['L_target_eff_'        m1 '_' m2])         = L_target_eff;
            DATA.(condName).meta.(['Ltot_'                m1 '_' m2 '_left']) = NtotalSide.(condName).left;
            DATA.(condName).meta.(['Ltot_'                m1 '_' m2 '_right'])= NtotalSide.(condName).right;
            DATA.(condName).meta.RNG.Type  = seedState.Type;
            DATA.(condName).meta.RNG.Seed  = seedState.Seed;

            % Seuil statistique et sommes brutes (pour recalcul ultérieur si besoin)
            DATA.(condName).(side_str).(['Seuil_' m1 '_' m2])  = seuil;
            DATA.(condName).(side_str).(['PxxSum_' m1 '_' m2]) = Pxx_sum;        % Σ|TFR1|²
            DATA.(condName).(side_str).(['PyySum_' m1 '_' m2]) = Pyy_sum;        % Σ|TFR2|²
            DATA.(condName).(side_str).(['PxySum_' m1 '_' m2]) = Pxy_sum;        % ΣTFr1·TFR2* (complexe)
            DATA.(condName).(side_str).(['L_' m1 '_' m2])      = L_cycles;
            DATA.(condName).(side_str).Freq = Freq;                               % Axe fréquentiel (Hz)

            % --- Génération des figures de contrôle qualité ---
            % 8 sous-figures : spectres de puissance, cohérence (significative et
            % profil fréquentiel) pour le cycle complet et les 4 sous-phases
            Time = linspace(0, 100, 1000); % Axe temporel normalisé (% du cycle)
            muscle1_name = m1;
            muscle2_name = m2;

            figure('Position',[100,100,1400,800]);

            % Spectre de puissance — Muscle 1
            subplot(2,4,1);
            imagesc(Time, Freq, PowSpec_s1);
            title(['Power Spec ' muscle1_name]); xlabel('% Cycle'); ylabel('Fréquence (Hz)'); colorbar;

            % Spectre de puissance — Muscle 2
            subplot(2,4,2);
            imagesc(Time, Freq, PowSpec_s2);
            title(['Power Spec ' muscle2_name]); xlabel('% Cycle'); ylabel('Fréquence (Hz)'); colorbar;

            % Cohérence significative (cohérence > seuil = 1, sinon 0)
            % avec repères visuels des événements de marche
            subplot(2,4,3);
            imagesc(Time, Freq, Coherence > seuil);
            title([muscle1_name '-' muscle2_name ' Coherence']); hold on;
            line([ge.opposite_toeoff     ge.opposite_toeoff],     [min(Freq) max(Freq)], 'Color','blue',  'LineWidth',2); % TO controlatéral
            line([ge.opposite_heelstrike ge.opposite_heelstrike], [min(Freq) max(Freq)], 'Color','green', 'LineWidth',2); % HS controlatéral
            line([ge.main_toeoff         ge.main_toeoff],         [min(Freq) max(Freq)], 'Color','red',   'LineWidth',2); % TO principal
            xlabel('% Cycle'); ylabel('Fréquence (Hz)');

            % Profil fréquentiel de la cohérence — Cycle complet
            subplot(2,4,4);
            plot(Freq, mean(Coherence,2)); hold on;
            line([0 400],[seuil seuil],'Color','red'); % Seuil de significativité
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Cycle complet'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            % Profil fréquentiel — Loading Response
            subplot(2,4,5);
            if ~isempty(loading_response_indices), plot(Freq, mean(Coherence(:,loading_response_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Loading Response'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            % Profil fréquentiel — Mid-Stance
            subplot(2,4,6);
            if ~isempty(midstance_indices), plot(Freq, mean(Coherence(:,midstance_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Mid Stance'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            % Profil fréquentiel — Pre-Swing
            subplot(2,4,7);
            if ~isempty(preswing_indices), plot(Freq, mean(Coherence(:,preswing_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Pre-swing'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            % Profil fréquentiel — Swing
            subplot(2,4,8);
            if ~isempty(swing_phase_indices), plot(Freq, mean(Coherence(:,swing_phase_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Swing'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            sgtitle(['Condition : ' condName ' | Côté : ' side_str ' | ' muscle1_name ' - ' muscle2_name]);

            % Sauvegarde de la figure au format PNG
            fig_name = sprintf('%s_%s_%s-%s', condName, side_str, muscle1_name, muscle2_name);
            saveas(gcf, fullfile(output_dir, [fig_name '.png']));
            close(gcf);

            % Affichage des résultats clés dans la console
            fprintf('\n=== Résultats pour %s-%s | %s | %s ===\n', muscle1_name, muscle2_name, condName, side_str);
            fprintf('  Opp TO : %.1f%% | Opp HS : %.1f%% | Main TO : %.1f%%\n', ...
                ge.opposite_toeoff, ge.opposite_heelstrike, ge.main_toeoff);

        end % Conditions
    end % Pairs

    % --- Sauvegarde finale du fichier de résultats ---
    save(fullfile(output_dir, ['Coherence_' Participant{iP} '.mat']), 'DATA');
    fprintf('Script terminé avec succès pour %s !\n', Participant{iP});

end % Participants