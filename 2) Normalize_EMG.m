%% NORMALIZE_EMG.m
%
% OBJECTIF
% --------
% Ce script extrait des cycles EMG (enveloppe nettoyée) depuis des fichiers .c3d,
% détecte des cycles aberrants (outliers), puis construit :
%   (1) une matrice de synergies "NNMF-ready" (SYNERGY_MATRIX)
%   (2) des cycles normalisés (CYCLES_SIGNAL_NORMALIZED)
%   (3) des cycles moyens (bruts et normalisés) et des infos associées (outliers, toe-off, nb cycles)
%
% La logique de normalisation est "in-task" (sans MVC), CEDE-compatible :
%   - On calcule une référence ref_max_plat *uniquement sur Plat* pour chaque
%     Participant × Jambe × Muscle :
%         ref_max_plat = max( tous les points de tous les cycles valides sur Plat )
%   - Cette référence est ensuite utilisée pour normaliser les cycles et la matrice de synergies
%     sur Plat, Medium et High.
%
% ENTRÉES
% -------
% - Fichiers .c3d nommés : <Participant>_<Condition>_<Essai>.c3d
%   ex : CTL_63_Plat_01.c3d
% - Données analogiques EMG dans les .c3d (via BTK)
% - Association capteur -> muscle :
%       Association.m  (doit définir sensor_association_left et sensor_association_right)
% - Fonctions de prétraitement EMG et détection événements :
%       filtrage, characterizeArtifacts, cleanEMGdata
%       indiceLeft / indiceRight, indiceLeftTO / indiceRightTO
% - Détection des outliers :
%       detect_outliers
%
% SORTIES
% -------
% Le script sauvegarde un .mat contenant notamment :
% - SYNERGY_MATRIX.(pid).(cond).(leg)                : [nCyclesValides*100 x nMuscles]
% - CYCLES_SIGNAL_NORMALIZED.(pid).(cond).(leg).(m)  : cycles normalisés (cell array)
% - CYCLES_SIGNAL.(pid).(cond).(leg).(m)             : cycles bruts (cell array)
% - CYCLES_MOYENS_BRUTS.(pid).Plat.(leg).(m)         : cycle moyen brut Plat
% - CYCLES_MOYENS.(pid).(cond).(leg).(m)             : cycle moyen normalisé (pour figures)
% - CYCLES_OUTLIERS.(pid).(cond).(leg).(m)           : outliers par muscle
% - CYCLES_COUNT.(pid).(cond).(leg)                  : nb cycles valides par condition/jambe
% - CYCLES_TOEOFF.(pid).(cond).(leg)                 : toe-off (% cycle) mean/std + distribution
% - REF_MAX_PLAT.(pid).(leg).(m)                     : référence de normalisation unique Plat
% - META_NORMALIZATION                               : métadonnées méthode
%
% NOTES MÉTHODOLOGIQUES IMPORTANTES
% --------------------------------
% - Les cycles sont définis HS -> HS (length(HS)-1 cycles)
% - Interpolation de chaque cycle sur 100 points (pchip)
% - Les outliers globaux d'une condition = UNION des cycles outliers sur tous les muscles
%   -> garantit que la SYNERGY_MATRIX aligne les mêmes cycles valides pour tous les muscles
% - Le toe-off % est calculé une seule fois par essai (m==1), pour éviter redondance

%% ============================ INITIALISATION ============================
clc; clear; close all;

% Ajoute le dossier du projet + sous-dossiers au path MATLAB
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire'));

% Se placer dans le répertoire contenant les .c3d à traiter
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Data\jeunes_enfants\');

%% ============================== PARAMÈTRES ==============================
Participant = {'CTL_63'};
Condition   = {'Plat','Medium','High'};
Essai       = {'01','02','03','04'}; % '05' '06' ...
muscles     = {'EMG_TAprox','EMG_TAdist','EMG_SOL','EMG_GM','EMG_VL','EMG_RF','EMG_ST','EMG_GMED'};
jambes      = {'left','right'};

% Mapping capteur->muscle (doit créer sensor_association_left/right)
run Association.m

%% =========================== STRUCTURES DE SORTIE ===========================
CYCLES_MOYENS_BRUTS      = struct();
CYCLES_MOYENS            = struct();
CYCLES_STD               = struct();
CYCLES_SIGNAL            = struct();
CYCLES_SIGNAL_NORMALIZED = struct();
CYCLES_OUTLIERS          = struct();
CYCLES_COUNT             = struct();
CYCLES_TOEOFF            = struct();
SYNERGY_MATRIX           = struct();

% Références de normalisation (unique par participant×jambe×muscle)
REF_MAX_PLAT             = struct();

%% ========================= MODE BATCH (FIGURES OFF) =========================
% Désactive l'affichage des figures pendant le traitement pour accélérer
% (les figures sont quand même sauvegardées sur disque)
original_visible = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

fprintf('=== DÉBUT DU TRAITEMENT ===\n');

%% =========================================================================
% BOUCLE PRINCIPALE
%   Jambe -> Participant -> Étape 1 (Plat pour la référence) -> Étape 2 (toutes conditions)
% =========================================================================
for j = 1:length(jambes)
    fprintf('Traitement jambe: %s\n', jambes{j});

    % ---------------------- ASSOCIATION CAPTEURS PAR JAMBE ----------------------
    % sensor_association.(muscle) doit retourner le nom du canal analogique correspondant
    if strcmp(jambes{j}, 'left')
        sensor_association = sensor_association_left;
    else
        sensor_association = sensor_association_right;
    end

    for iP = 1:length(Participant)
        pid = Participant{iP};
        fprintf('  Participant: %s\n', pid);

        % ===================== ÉTAPE 1 — PLAT (RÉFÉRENCE UNIQUE) =====================
        % Objectif : extraire cycles bruts Plat + détecter outliers + définir ref_max_plat
        fprintf('  === ÉTAPE 1: Traitement de PLAT (référence unique) ===\n');

        iC_plat = find(strcmp(Condition,'Plat'),1);
        if isempty(iC_plat)
            error('Condition "Plat" introuvable dans la liste Condition.');
        end
        cond_ref = Condition{iC_plat};

        % Dossier de sortie des figures "cycles concaténés"
        participant_folder_plat = fullfile( ...
            'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ...
            pid, cond_ref);
        if ~exist(participant_folder_plat,'dir'); mkdir(participant_folder_plat); end

        fprintf('    Condition: %s (référence)\n', cond_ref);
        fprintf('      Extraction des cycles Plat + Toe-Off...\n');

        % On remplit :
        %   CYCLES_SIGNAL.(pid).Plat.(leg).(muscle) = {cycles bruts interpolés (100pts)}
        %   CYCLES_OUTLIERS.(...)                   = outliers (par muscle)
        %   CYCLES_TOEOFF.(...)                     = stats toe-off (% cycle)
        for m = 1:length(muscles)
            muscle = muscles{m};

            % Figure de contrôle : tous les cycles concaténés dans un même graphe
            figure_concatenated = figure('Name', sprintf('Cycles Concaténés - %s - %s - %s - %s', pid, cond_ref, jambes{j}, muscle));
            hold on;

            % Accumulateurs pour 1 muscle sur l'ensemble des essais
            cycle_offset    = 0;     % index "temps concaténé" pour tracer les cycles bout à bout
            cycle_positions = [];    % positions où un cycle se termine (pour les traits rouges)
            muscle_maxima   = [];    % max par cycle (sert à detect_outliers)
            total_cycles    = 0;     % compteur de cycles extraits
            all_cycles_data = {};    % cycles bruts interpolés (100 points) stockés en cellule

            % Toe-off (% cycle) : calculé uniquement quand m==1 (évite redondance)
            all_toeoff_percentages = [];

            for iEs = 1:length(Essai)
                file = [pid '_' cond_ref '_' Essai{iEs} '.c3d'];

                try
                    % Lecture .c3d et analogiques
                    data    = btkReadAcquisition(file);
                    analogs = btkGetAnalogs(data);

                    % Sélection du canal EMG associé au muscle
                    sensor_name = sensor_association.(muscle);
                    EMG_signal  = analogs.(sensor_name);

                    % -------------------------- PRÉTRAITEMENT EMG --------------------------
                    % FreqS : fréquence analogique (EMG)
                    % FreqVicon : fréquence events (ici fixée à 100 Hz)
                    FreqS      = btkGetAnalogFrequency(data);
                    FreqVicon  = 100;

                    % Filtrage bande passante EMG + nettoyage + enveloppe
                    % (les paramètres exacts dépendent de tes fonctions)
                    EMGfilt = filtrage(EMG_signal, FreqS, 20, 400);
                    artifacts_info = characterizeArtifacts(EMGfilt, FreqS); 
                    [signal_cleaned, EMGenvcleaned] = cleanEMGdata(EMGfilt, FreqS, 'STD', 5, 15, false);

                    % -------------------------- ÉVÉNEMENTS GAIT --------------------------
                    % HS / TO indexés en "frames" compatibles avec l'EMG (selon tes fonctions indice*)
                    if strcmp(jambes{j}, 'left')
                        HS = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_signal);
                        TO = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                    else
                        HS = indiceRight(data, analogs, FreqS, FreqVicon, EMG_signal);
                        TO = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                    end

                    % Définition cycles : HS(i) -> HS(i+1)
                    num_cycles = length(HS)-1;
                    cycles_frames = cell(num_cycles,1);

                    % -------------------------- TOE-OFF % CYCLE --------------------------
                    % Calculé 1 seule fois (m==1) car indépendant du muscle
                    if m == 1
                        for ii = 1:num_cycles
                            cstart = HS(ii); cend = HS(ii+1);
                            dur = cend - cstart;

                            TO_in = TO(TO > cstart & TO < cend);
                            if ~isempty(TO_in)
                                TO_frame = TO_in(1);
                                TO_percentage = ((TO_frame - cstart)/dur)*100;
                                all_toeoff_percentages = [all_toeoff_percentages, TO_percentage];
                            else
                                % Si TO absent : tu mets 60% par défaut (hypothèse “toe-off ~60%”)
                                warning('Pas de toe-off (Plat) : cycle %d, essai %s', ii, Essai{iEs});
                                all_toeoff_percentages = [all_toeoff_percentages, 60];
                            end
                        end
                    end

                    % Construire les indices (frames) de chaque cycle
                    for ii = 1:num_cycles
                        cycles_frames{ii} = HS(ii):HS(ii+1);
                    end

                    % --------------------- EXTRACTION + INTERPOLATION (100 pts) ---------------------
                    for ii = 1:num_cycles
                        cdata   = EMGenvcleaned(cycles_frames{ii});
                        cinterp = interp1(linspace(1,length(cdata),length(cdata)), cdata, ...
                                          linspace(1,length(cdata),100), 'pchip');

                        total_cycles = total_cycles + 1;
                        all_cycles_data{total_cycles} = cinterp;

                        % Pour QA/visualisation : cycles concaténés
                        cycle_positions = [cycle_positions, cycle_offset + length(cinterp)];
                        muscle_maxima   = [muscle_maxima,   max(cinterp)];
                        plot(cycle_offset + (1:length(cinterp)), cinterp, 'b');

                        cycle_offset = cycle_offset + length(cinterp);
                    end

                catch ME
                    warning('Erreur fichier %s: %s', file, ME.message);
                    continue;
                end
            end

            % Séparateurs entre cycles (traits verticaux rouges)
            y_limits = ylim;
            for ii = 1:length(cycle_positions)
                plot([cycle_positions(ii), cycle_positions(ii)], y_limits, 'r', 'LineWidth', 0.5);
            end
            title(sprintf('All cycles for %s - %s - Muscle: %s - Leg: %s', pid, cond_ref, muscle, jambes{j}), 'Interpreter','none');
            xlabel('Temps normalisé'); ylabel('Enveloppe');

            % -------------------------- OUTLIERS (PAR MUSCLE) --------------------------
            % detect_outliers : attend un indicateur par cycle (ici max(cycle))
            % -> retourne typiquement un vecteur logique (1 = outlier)
            outlier_indices = detect_outliers(muscle_maxima, cycle_offset, 3);

            % Sauvegarde de la figure QA
            save_path = fullfile(participant_folder_plat, sprintf('%s_%s_%s_%s.png', pid, cond_ref, jambes{j}, muscle));
            print(figure_concatenated, save_path, '-dpng', '-r150'); close(figure_concatenated);

            % -------------------------- STOCKAGE CYCLES + OUTLIERS --------------------------
            if ~isfield(CYCLES_SIGNAL, pid); CYCLES_SIGNAL.(pid)=struct(); end
            if ~isfield(CYCLES_SIGNAL.(pid), cond_ref); CYCLES_SIGNAL.(pid).(cond_ref)=struct(); end
            if ~isfield(CYCLES_SIGNAL.(pid).(cond_ref), jambes{j}); CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j})=struct(); end
            CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j}).(muscle) = all_cycles_data;

            if ~isfield(CYCLES_OUTLIERS, pid); CYCLES_OUTLIERS.(pid)=struct(); end
            if ~isfield(CYCLES_OUTLIERS.(pid), cond_ref); CYCLES_OUTLIERS.(pid).(cond_ref)=struct(); end
            if ~isfield(CYCLES_OUTLIERS.(pid).(cond_ref), jambes{j}); CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j})=struct(); end
            CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j}).(muscle) = outlier_indices;

            % -------------------------- STOCKAGE TOE-OFF (% CYCLE) -------------------------
            if ~isempty(all_toeoff_percentages)
                if ~isfield(CYCLES_TOEOFF, pid); CYCLES_TOEOFF.(pid)=struct(); end
                if ~isfield(CYCLES_TOEOFF.(pid), cond_ref); CYCLES_TOEOFF.(pid).(cond_ref)=struct(); end
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).mean_percentage = mean(all_toeoff_percentages);
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).std_percentage  = std(all_toeoff_percentages);
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).all_percentages = all_toeoff_percentages;
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).num_cycles      = length(all_toeoff_percentages);
            end
        end

        % -------------------------- OUTLIERS GLOBAUX PLAT (UNION MUSCLES) --------------------------
        fprintf('      Identification des outliers globaux Plat...\n');
        all_outlier_indices_plat = [];
        for m = 1:length(muscles)
            muscle = muscles{m};
            if isfield(CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j}), muscle)
                outlier_idx = CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j}).(muscle);
                all_outlier_indices_plat = union(all_outlier_indices_plat, find(outlier_idx));
            end
        end

        % -------------------------- CYCLE MOYEN BRUT + REF_MAX_PLAT --------------------------
        fprintf('      Calcul des cycles moyens BRUTS + REF_MAX_PLAT (max sur tous cycles valides)...\n');
        for m = 1:length(muscles)
            muscle = muscles{m};
            if ~isfield(CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j}), muscle), continue; end
            all_cycles = CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j}).(muscle);

            % Cycles valides (on exclut l’union des outliers)
            if ~isempty(all_outlier_indices_plat)
                total_avail = length(all_cycles);
                valid_idx   = setdiff(1:total_avail, all_outlier_indices_plat);
            else
                valid_idx   = 1:length(all_cycles);
            end
            if isempty(valid_idx), continue; end

            % Matrice [nCyclesValides x 100] pour faire mean et max global
            all_points = zeros(length(valid_idx),100);
            for k = 1:length(valid_idx)
                cnum = valid_idx(k);
                if cnum <= length(all_cycles)
                    all_points(k,:) = all_cycles{cnum};
                end
            end

            cycle_moyen_brut = mean(all_points,1);

            % Référence unique : max global sur tous points de tous cycles valides (Plat)
            ref_max_plat = max(all_points(:));
            if ref_max_plat <= 0
                warning('ref_max_plat <= 0 pour %s - %s - %s. Fallback: max du cycle moyen.', pid, jambes{j}, muscle);
                ref_max_plat = max(cycle_moyen_brut);
            end

            % Stockages
            if ~isfield(CYCLES_MOYENS_BRUTS, pid), CYCLES_MOYENS_BRUTS.(pid)=struct(); end
            if ~isfield(CYCLES_MOYENS_BRUTS.(pid), cond_ref), CYCLES_MOYENS_BRUTS.(pid).(cond_ref)=struct(); end
            if ~isfield(CYCLES_MOYENS_BRUTS.(pid).(cond_ref), jambes{j}), CYCLES_MOYENS_BRUTS.(pid).(cond_ref).(jambes{j})=struct(); end
            CYCLES_MOYENS_BRUTS.(pid).(cond_ref).(jambes{j}).(muscle) = cycle_moyen_brut;

            if ~isfield(REF_MAX_PLAT, pid), REF_MAX_PLAT.(pid)=struct(); end
            if ~isfield(REF_MAX_PLAT.(pid), jambes{j}), REF_MAX_PLAT.(pid).(jambes{j})=struct(); end
            REF_MAX_PLAT.(pid).(jambes{j}).(muscle) = ref_max_plat;

            fprintf('        %s: REF_MAX_PLAT = %.6f\n', muscle, ref_max_plat);
        end

        % ===================== ÉTAPE 2 — TOUTES CONDITIONS =====================
        % Objectif : pour chaque condition :
        %   - extraire cycles bruts (si cond ≠ Plat)
        %   - calculer outliers globaux (union muscles)
        %   - construire SYNERGY_MATRIX normalisée par ref_max_plat
        %   - exporter cycles normalisés + cycle moyen/SD (figures)
        fprintf('  === ÉTAPE 2: Construction des matrices avec normalisation unique ===\n');

        for iC = 1:length(Condition)
            cond = Condition{iC};
            fprintf('    Condition: %s\n', cond);

            participant_folder = fullfile( ...
                'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ...
                pid, cond);
            if ~exist(participant_folder,'dir'); mkdir(participant_folder); end

            % -------------------------- EXTRACTION CYCLES (si cond ≠ Plat) --------------------------
            % Pour Plat, les cycles ont déjà été extraits à l’étape 1.
            if ~strcmp(cond,'Plat')
                fprintf('      Extraction cycles + toe-off (%s)...\n', cond);
                for m = 1:length(muscles)
                    muscle = muscles{m};

                    figure_concatenated = figure('Name', sprintf('Cycles Concaténés - %s - %s - %s - %s', pid, cond, jambes{j}, muscle));
                    hold on;

                    cycle_offset    = 0;
                    cycle_positions = [];
                    muscle_maxima   = [];
                    total_cycles    = 0;
                    all_cycles_data = {};
                    all_toeoff_percentages = [];

                    for iEs = 1:length(Essai)
                        file = [pid '_' cond '_' Essai{iEs} '.c3d'];
                        try
                            data    = btkReadAcquisition(file);
                            analogs = btkGetAnalogs(data);

                            sensor_name = sensor_association.(muscle);
                            EMG_signal  = analogs.(sensor_name);

                            FreqS      = btkGetAnalogFrequency(data);
                            FreqVicon  = 100;

                            EMGfilt = filtrage(EMG_signal, FreqS, 20, 400);
                            artifacts_info = characterizeArtifacts(EMGfilt, FreqS);
                            [signal_cleaned, EMGenvcleaned] = cleanEMGdata(EMGfilt, FreqS, 'STD', 5, 15, false);

                            if strcmp(jambes{j}, 'left')
                                HS = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_signal);
                                TO = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                            else
                                HS = indiceRight(data, analogs, FreqS, FreqVicon, EMG_signal);
                                TO = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                            end

                            num_cycles = length(HS)-1;
                            cycles_frames = cell(num_cycles,1);

                            % Toe-off % cycle (m==1 uniquement)
                            if m == 1
                                for ii = 1:num_cycles
                                    cstart = HS(ii); cend = HS(ii+1);
                                    dur = cend - cstart;
                                    TO_in = TO(TO > cstart & TO < cend);
                                    if ~isempty(TO_in)
                                        TO_frame = TO_in(1);
                                        TO_percentage = ((TO_frame - cstart)/dur)*100;
                                        all_toeoff_percentages = [all_toeoff_percentages, TO_percentage];
                                    else
                                        warning('Pas de toe-off (%s) : cycle %d, essai %s', cond, ii, Essai{iEs});
                                        all_toeoff_percentages = [all_toeoff_percentages, 60];
                                    end
                                end
                            end

                            for ii = 1:num_cycles
                                cycles_frames{ii} = HS(ii):HS(ii+1);
                            end

                            for ii = 1:num_cycles
                                cdata   = EMGenvcleaned(cycles_frames{ii});
                                cinterp = interp1(linspace(1,length(cdata),length(cdata)), cdata, ...
                                                  linspace(1,length(cdata),100), 'pchip');

                                total_cycles = total_cycles + 1;
                                all_cycles_data{total_cycles} = cinterp;

                                cycle_positions = [cycle_positions, cycle_offset + length(cinterp)];
                                muscle_maxima   = [muscle_maxima,   max(cinterp)];

                                plot(cycle_offset + (1:length(cinterp)), cinterp, 'b');
                                cycle_offset = cycle_offset + length(cinterp);
                            end

                        catch ME
                            warning('Erreur fichier %s: %s', file, ME.message);
                            continue;
                        end
                    end

                    y_limits = ylim;
                    for ii = 1:length(cycle_positions)
                        plot([cycle_positions(ii), cycle_positions(ii)], y_limits, 'r', 'LineWidth', 0.5);
                    end
                    title(sprintf('All cycles for %s - %s - Muscle: %s - Leg: %s', pid, cond, muscle, jambes{j}), 'Interpreter','none');
                    xlabel('Temps normalisé'); ylabel('Enveloppe');

                    outlier_indices = detect_outliers(muscle_maxima, cycle_offset, 3);

                    save_path = fullfile(participant_folder, sprintf('%s_%s_%s_%s.png', pid, cond, jambes{j}, muscle));
                    print(figure_concatenated, save_path, '-dpng', '-r150'); close(figure_concatenated);

                    if ~isfield(CYCLES_SIGNAL, pid), CYCLES_SIGNAL.(pid)=struct(); end
                    if ~isfield(CYCLES_SIGNAL.(pid), cond), CYCLES_SIGNAL.(pid).(cond)=struct(); end
                    if ~isfield(CYCLES_SIGNAL.(pid).(cond), jambes{j}), CYCLES_SIGNAL.(pid).(cond).(jambes{j})=struct(); end
                    CYCLES_SIGNAL.(pid).(cond).(jambes{j}).(muscle) = all_cycles_data;

                    if ~isfield(CYCLES_OUTLIERS, pid), CYCLES_OUTLIERS.(pid)=struct(); end
                    if ~isfield(CYCLES_OUTLIERS.(pid), cond), CYCLES_OUTLIERS.(pid).(cond)=struct(); end
                    if ~isfield(CYCLES_OUTLIERS.(pid).(cond), jambes{j}), CYCLES_OUTLIERS.(pid).(cond).(jambes{j})=struct(); end
                    CYCLES_OUTLIERS.(pid).(cond).(jambes{j}).(muscle) = outlier_indices;

                    if ~isempty(all_toeoff_percentages)
                        if ~isfield(CYCLES_TOEOFF, pid), CYCLES_TOEOFF.(pid)=struct(); end
                        if ~isfield(CYCLES_TOEOFF.(pid), cond), CYCLES_TOEOFF.(pid).(cond)=struct(); end
                        CYCLES_TOEOFF.(pid).(cond).(jambes{j}).mean_percentage = mean(all_toeoff_percentages);
                        CYCLES_TOEOFF.(pid).(cond).(jambes{j}).std_percentage  = std(all_toeoff_percentages);
                        CYCLES_TOEOFF.(pid).(cond).(jambes{j}).all_percentages = all_toeoff_percentages;
                        CYCLES_TOEOFF.(pid).(cond).(jambes{j}).num_cycles      = length(all_toeoff_percentages);
                    end
                end
            end

            % -------------------------- OUTLIERS GLOBAUX (UNION MUSCLES) --------------------------
            all_outlier_indices = [];
            for m = 1:length(muscles)
                muscle = muscles{m};
                if isfield(CYCLES_OUTLIERS.(pid).(cond).(jambes{j}), muscle)
                    out_idx = CYCLES_OUTLIERS.(pid).(cond).(jambes{j}).(muscle);
                    all_outlier_indices = union(all_outlier_indices, find(out_idx));
                end
            end

            % Déterminer nb cycles total (basé sur le 1er muscle disponible)
            if isfield(CYCLES_SIGNAL.(pid).(cond).(jambes{j}), muscles{1})
                total_avail = length(CYCLES_SIGNAL.(pid).(cond).(jambes{j}).(muscles{1}));
            else
                total_avail = 0;
            end

            % Indices cycles valides = tous - outliers globaux
            if total_avail>0 && ~isempty(all_outlier_indices)
                valid_idx = setdiff(1:total_avail, all_outlier_indices);
            elseif total_avail>0
                valid_idx = 1:total_avail;
            else
                valid_idx = [];
            end

            num_valid_cycles = length(valid_idx);

            if ~isfield(CYCLES_COUNT, pid), CYCLES_COUNT.(pid)=struct(); end
            if ~isfield(CYCLES_COUNT.(pid), cond), CYCLES_COUNT.(pid).(cond)=struct(); end
            CYCLES_COUNT.(pid).(cond).(jambes{j}) = num_valid_cycles;

            fprintf('        Nombre de cycles valides pour %s: %d\n', cond, num_valid_cycles);
            if num_valid_cycles == 0
                warning('Aucun cycle valide: %s - %s - %s', pid, cond, jambes{j});
                continue;
            end

            % -------------------------- SYNERGY_MATRIX (NNMF-ready) --------------------------
            % Matrice finale :
            %   lignes = nCyclesValides * 100 points
            %   colonnes = nMuscles
            % Chaque colonne correspond à la concaténation des cycles normalisés (100 pts) du muscle.
            synergy_matrix = zeros(num_valid_cycles*100, length(muscles));

            for m = 1:length(muscles)
                muscle = muscles{m};

                if ~isfield(CYCLES_SIGNAL.(pid).(cond).(jambes{j}), muscle), continue; end
                all_cycles = CYCLES_SIGNAL.(pid).(cond).(jambes{j}).(muscle);

                % Référence unique (Plat) pour cette jambe+muscle
                if isfield(REF_MAX_PLAT,pid) && isfield(REF_MAX_PLAT.(pid),jambes{j}) && isfield(REF_MAX_PLAT.(pid).(jambes{j}),muscle)
                    ref_max_plat = REF_MAX_PLAT.(pid).(jambes{j}).(muscle);
                else
                    warning('REF_MAX_PLAT manquant pour %s - %s - %s. Fallback: max individuel.', pid, jambes{j}, muscle);
                    ref_max_plat = NaN;
                end

                % Construire la colonne du muscle : concat cycles normalisés
                muscle_column = [];
                for kk = 1:length(valid_idx)
                    cnum = valid_idx(kk);
                    if cnum <= length(all_cycles)
                        cdata = all_cycles{cnum}; % 100 points
                        if ~isnan(ref_max_plat) && ref_max_plat > 0
                            cdata_norm = cdata / ref_max_plat; % normalisation unique
                        else
                            % fallback : normalisation cycle par cycle
                            cmax = max(cdata);
                            if cmax>0
                                cdata_norm = cdata / cmax;
                            else
                                cdata_norm = cdata;
                            end
                        end
                        muscle_column = [muscle_column; cdata_norm(:)];
                    end
                end

                % Remplissage matrice
                if length(muscle_column) == num_valid_cycles*100
                    synergy_matrix(:,m) = muscle_column;
                else
                    warning('Taille incorrecte (%s - %s - %s): attendu %d, obtenu %d', muscle, pid, cond, num_valid_cycles*100, length(muscle_column));
                end

                % -------------------------- CYCLE MOYEN + SD (NORMALISÉ) --------------------------
                % Utilisé uniquement pour générer des figures de synthèse par muscle/cond/jambe
                n_cycles = num_valid_cycles; n_points = 100;
                if length(muscle_column) == n_cycles*n_points
                    muscle_matrix = reshape(muscle_column,[n_points,n_cycles])';
                    mean_cycle = mean(muscle_matrix,1);
                    sd_cycle   = std(muscle_matrix,0,1);

                    if ~isfield(CYCLES_MOYENS, pid), CYCLES_MOYENS.(pid)=struct(); end
                    if ~isfield(CYCLES_MOYENS.(pid), cond), CYCLES_MOYENS.(pid).(cond)=struct(); end
                    if ~isfield(CYCLES_MOYENS.(pid).(cond), jambes{j}), CYCLES_MOYENS.(pid).(cond).(jambes{j})=struct(); end
                    CYCLES_MOYENS.(pid).(cond).(jambes{j}).(muscle) = mean_cycle;

                    % Figure sauvegardée (mean ± SD)
                    fig = figure('Visible','off'); hold on;
                    x = linspace(0,100,n_points);
                    fill([x,fliplr(x)], [mean_cycle+sd_cycle, fliplr(mean_cycle - sd_cycle)], ...
                         [1 0 0], 'FaceAlpha',0.3,'EdgeColor','none');
                    plot(x, mean_cycle, 'r','LineWidth',2);
                    title(sprintf('Cycle moyen normalisé - %s - %s - %s - %s', pid, cond, jambes{j}, muscle), 'Interpreter','none');
                    xlabel('% du cycle'); ylabel('Activation normalisée'); xlim([0 100]);

                    fig_folder = fullfile(participant_folder,'CycleMoyen');
                    if ~exist(fig_folder,'dir'); mkdir(fig_folder); end
                    fig_path = fullfile(fig_folder, sprintf('%s_%s_%s_%s_cycleMoyen.png', pid, cond, jambes{j}, muscle));
                    print(fig, fig_path, '-dpng', '-r150'); close(fig);
                end

                % -------------------------- EXPORT CYCLES NORMALISÉS --------------------------
                normalized_cycles = cell(1,length(valid_idx));
                for kk = 1:length(valid_idx)
                    cnum = valid_idx(kk);
                    if cnum <= length(all_cycles)
                        cdata = all_cycles{cnum};
                        if ~isnan(ref_max_plat) && ref_max_plat > 0
                            normalized_cycles{kk} = cdata / ref_max_plat;
                        else
                            cmax = max(cdata);
                            if cmax > 0
                                normalized_cycles{kk} = cdata / cmax;
                            else
                                normalized_cycles{kk} = cdata;
                            end
                        end
                    end
                end

                if ~isfield(CYCLES_SIGNAL_NORMALIZED, pid), CYCLES_SIGNAL_NORMALIZED.(pid)=struct(); end
                if ~isfield(CYCLES_SIGNAL_NORMALIZED.(pid), cond), CYCLES_SIGNAL_NORMALIZED.(pid).(cond)=struct(); end
                if ~isfield(CYCLES_SIGNAL_NORMALIZED.(pid).(cond), jambes{j}), CYCLES_SIGNAL_NORMALIZED.(pid).(cond).(jambes{j})=struct(); end
                CYCLES_SIGNAL_NORMALIZED.(pid).(cond).(jambes{j}).(muscle) = normalized_cycles;

            end % muscles

            % Stockage matrice synergies
            if ~isfield(SYNERGY_MATRIX, pid), SYNERGY_MATRIX.(pid)=struct(); end
            if ~isfield(SYNERGY_MATRIX.(pid), cond), SYNERGY_MATRIX.(pid).(cond)=struct(); end
            SYNERGY_MATRIX.(pid).(cond).(jambes{j}) = synergy_matrix;

            fprintf('        SYNERGY_MATRIX: %d x %d (points concaténés x muscles)\n', size(synergy_matrix,1), size(synergy_matrix,2));
        end % Condition
    end % Participant
end % jambes

% ============================ RESTORE FIGURES ============================
set(0,'DefaultFigureVisible', original_visible);

%% =============================== SAUVEGARDE ===============================
fprintf('\n=== SAUVEGARDE DES RÉSULTATS ===\n');

save_dir  = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS';
if ~exist(save_dir,'dir'); mkdir(save_dir); end

% NB : Participant{iP} ici correspond au dernier iP de la boucle (ok si un seul participant).
% Si plusieurs participants, préférer sauvegarder par pid dans la boucle, ou concaténer.
save_file = fullfile(save_dir, [Participant{iP} '_MATRIX.mat']);

% Métadonnées : trace la méthode de normalisation
META_NORMALIZATION.method               = 'Single in-task reference (Plat)';
META_NORMALIZATION.reference_definition = 'ref_max_plat = max over all points of all valid cycles in Plat (per participant×leg×muscle)';
META_NORMALIZATION.applied_to           = 'SYNERGY_MATRIX and CYCLES_SIGNAL_NORMALIZED for Plat/Medium/High';
META_NORMALIZATION.notes                = 'Outliers removed before computing ref; fallback per-cycle max used only if ref missing or <=0';

save(save_file, 'SYNERGY_MATRIX', 'CYCLES_SIGNAL_NORMALIZED', ...
     'CYCLES_MOYENS_BRUTS', 'CYCLES_MOYENS', 'CYCLES_COUNT', ...
     'CYCLES_OUTLIERS', 'CYCLES_TOEOFF', 'REF_MAX_PLAT', 'META_NORMALIZATION', '-v7.3');

fprintf('Structures sauvegardées dans: %s\n', save_file);
fprintf('\n=== TRAITEMENT TERMINÉ ===\n');