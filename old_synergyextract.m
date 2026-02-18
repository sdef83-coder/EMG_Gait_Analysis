%% EXTRAIRE LA MATRICE DE SYNERGY ET LES MATRICES INDIVIDUELLES
% Les matrices indiv. contiennent : 
% 1) Matrices de synergies; 
% 2) Signaux de chaque cycles normalisés; 
% 3) Cycles moyens bruts & normalisés/muscle/jambe/condition;
% 4) Les cycles outliers
% 5) Le nombre de cycle/condition/jambe
% 6) Les % du cycle des Toe-offs

clc;
clear;
close all;

addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire'));

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Data\\jeunes_enfants\');

Participant = {'CTL_63'};
Condition = {'Plat' 'Medium' 'High'};
Essai = {'01' '02' '03' '04'}; %'05' '06' '07' '08' '09' '10'
muscles = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'}; %'EMG_GMED' % Penser à retirer muscle si pas traité
jambes = {'left' 'right'};

run Association.m

% Initialisation des structures
CYCLES_MOYENS_BRUTS = struct();
CYCLES_MOYENS = struct();
CYCLES_STD = struct();
CYCLES_SIGNAL = struct();
CYCLES_SIGNAL_NORMALIZED = struct();
CYCLES_OUTLIERS = struct();
MEAN_MAX_PLAT = struct(); 
MUSCLE_DATA = struct();
SYNERGY_MATRIX = struct();
CYCLES_COUNT = struct();
MAX_CYCLES_BRUTS = struct();
CYCLES_TOEOFF = struct(); 

% Désactiver l'affichage des figures pendant le traitement
original_visible = get(0, 'DefaultFigureVisible');
set(0, 'DefaultFigureVisible', 'off');

fprintf('=== DÉBUT DU TRAITEMENT ===\n');

% Boucle des Jambes
for j = 1:length(jambes)
    fprintf('Traitement jambe: %s\n', jambes{j});
    
    % Association de capteurs
    if strcmp(jambes{j}, 'left')
        sensor_association = sensor_association_left;
    else
        sensor_association = sensor_association_right;
    end
    
    for iP = 1:length(Participant)
        fprintf('  Participant: %s\n', Participant{iP});
        
        % ÉTAPE PRÉLIMINAIRE: Traiter d'abord la condition Plat pour établir les références
        fprintf('  === ÉTAPE 1: Traitement de la condition Plat (référence) ===\n');
        
        for iC = 1:length(Condition)
            if strcmp(Condition{iC}, 'Plat')
                fprintf('    Condition: %s (référence)\n', Condition{iC});
                
                % Dossier d'enregistrement
                participant_folder = fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle\', Participant{iP}, Condition{iC});
                if ~exist(participant_folder, 'dir')
                    mkdir(participant_folder);
                end

                % Extraction et stockage de tous les cycles pour Plat
                fprintf('      Extraction des cycles Plat...\n');
                fprintf('      Extraction des toe-offs Plat...\n');
                all_toeoff_percentages = [];
                
                % Traitement muscle par muscle
                for m = 1:length(muscles)
                    muscle = muscles{m};
                    
                    % Créer figure pour visualisation
                    figure_concatenated = figure('Name', sprintf('Cycles Concaténés - %s - %s - %s - %s', Participant{iP}, Condition{iC}, jambes{j}, muscle));
                    hold on;
                    
                    cycle_offset = 0;
                    cycle_positions = [];
                    muscle_maxima = [];
                    total_cycles = 0;
                    all_cycles_data = {}; % Cell array pour stocker tous les cycles

                    % Boucle des essais
                    for iEs = 1:length(Essai)  
                        file = [Participant{iP} '_' Condition{iC} '_' Essai{iEs} '.c3d'];

                        try
                            % EXTRACTION DES DONNEES BRUTES
                            data = btkReadAcquisition(file);
                            analogs = btkGetAnalogs(data);
                            sensor_name = sensor_association.(muscle);
                            EMG_signal = analogs.(sensor_name);

                            % Traitement du signal
                            FreqS = btkGetAnalogFrequency(data);
                            FreqVicon = 100;
                            EMGfilt = filtrage(EMG_signal, FreqS, 20, 450);
                            artifacts_info = characterizeArtifacts(EMGfilt,FreqS);
                            [signal_cleaned, EMGenvcleaned] = cleanEMGdata(EMGfilt, FreqS, 'STD', 5, 15, false);

                            % Détection HS
                            if strcmp(jambes{j}, 'left')
                                HS = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_signal);
                                TO = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                            else
                                HS = indiceRight(data, analogs, FreqS, FreqVicon, EMG_signal);
                                TO = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                            end
                            
                            num_cycles = length(HS) - 1;
                            cycles_frames = cell(num_cycles, 1);

                            % === CALCUL DES POURCENTAGES TOE-OFF POUR CET ESSAI ===
                            if m == 1 % Ne calculer qu'une seule fois par essai (pour le premier muscle)
                                for i = 1:num_cycles
                                    cycle_start = HS(i);
                                    cycle_end = HS(i + 1);
                                    cycle_duration = cycle_end - cycle_start;
                                    
                                    % Trouver le TO dans ce cycle
                                    TO_in_cycle = TO(TO > cycle_start & TO < cycle_end);
                                    
                                    if ~isempty(TO_in_cycle)
                                        % Prendre le premier TO dans le cycle
                                        TO_frame = TO_in_cycle(1);
                                        % Calculer le pourcentage
                                        TO_percentage = ((TO_frame - cycle_start) / cycle_duration) * 100;
                                        all_toeoff_percentages = [all_toeoff_percentages, TO_percentage];
                                    else
                                        % Si pas de TO trouvé, utiliser une valeur par défaut (environ 60% du cycle)
                                        warning('Pas de toe-off trouvé pour le cycle %d, essai %s, condition %s', i, Essai{iEs}, Condition{iC});
                                        all_toeoff_percentages = [all_toeoff_percentages, 60]; % Valeur typique
                                    end
                                end
                            end

                            for i = 1:num_cycles
                                cycles_frames{i} = HS(i):HS(i + 1);
                            end

                            % Traitement des cycles individuels
                            for i = 1:num_cycles
                                cycle_signal = EMGenvcleaned(cycles_frames{i});
                                cycle_signal_interpolated = interp1(linspace(1, length(cycle_signal), length(cycle_signal)), cycle_signal, linspace(1, length(cycle_signal), 100), 'pchip');
                                
                                % Stocker le cycle BRUT (non normalisé)
                                total_cycles = total_cycles + 1;
                                all_cycles_data{total_cycles} = cycle_signal_interpolated;
                                
                                cycle_positions = [cycle_positions, cycle_offset + length(cycle_signal_interpolated)];
                                muscle_maxima = [muscle_maxima, max(cycle_signal_interpolated)];

                                % Tracé individuel de chaque cycle
                                plot(cycle_offset + (1:length(cycle_signal_interpolated)), cycle_signal_interpolated, 'b');
                                cycle_offset = cycle_offset + length(cycle_signal_interpolated);
                            end
                            
                        catch ME
                            warning('Erreur lors du traitement du fichier %s: %s', file, ME.message);
                            continue;
                        end
                    end

                    % Finaliser la figure
                    y_limits = ylim;
                    for i = 1:length(cycle_positions)
                        plot([cycle_positions(i), cycle_positions(i)], y_limits, 'r', 'LineWidth', 0.5);
                    end

                    title_str = sprintf('All cycles for %s - %s - Muscle: %s - Leg: %s', ...
                                        Participant{iP}, Condition{iC}, muscle, jambes{j});
                    title(title_str, 'Interpreter', 'none');
                    xlabel('Temps normalisé');
                    ylabel('Enveloppe');

                    % Détection outliers basée sur les maxima
                    outlier_indices = detect_outliers(muscle_maxima, cycle_offset, 3);

                    % Sauvegarde figure
                    save_path = fullfile(participant_folder, sprintf('%s_%s_%s_%s.png', Participant{iP}, Condition{iC}, jambes{j}, muscle));
                    print(figure_concatenated, save_path, '-dpng', '-r150');
                    close(figure_concatenated);

                    % Stockage des données et outliers
                    if ~isfield(CYCLES_SIGNAL, Participant{iP})
                        CYCLES_SIGNAL.(Participant{iP}) = struct();
                    end
                    if ~isfield(CYCLES_SIGNAL.(Participant{iP}), Condition{iC})
                        CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}) = struct();
                    end
                    if ~isfield(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}), jambes{j})
                        CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
                    end
                    
                    CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = all_cycles_data;
                    
                    % Structures outliers
                    if ~isfield(CYCLES_OUTLIERS, Participant{iP})
                        CYCLES_OUTLIERS.(Participant{iP}) = struct();
                    end
                    if ~isfield(CYCLES_OUTLIERS.(Participant{iP}), Condition{iC})
                        CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}) = struct();
                    end
                    if ~isfield(CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}), jambes{j})
                        CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
                    end
                    
                    CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = outlier_indices;
                end

                % === STOCKAGE TOE-OFF MOYEN POUR PLAT ===
                if ~isempty(all_toeoff_percentages)
                    mean_toeoff = mean(all_toeoff_percentages);
                    std_toeoff = std(all_toeoff_percentages);
                    
                    % Créer la structure CYCLES_TOEOFF
                    if ~isfield(CYCLES_TOEOFF, Participant{iP})
                        CYCLES_TOEOFF.(Participant{iP}) = struct();
                    end
                    if ~isfield(CYCLES_TOEOFF.(Participant{iP}), Condition{iC})
                        CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}) = struct();
                    end
                    
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).mean_percentage = mean_toeoff;
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).std_percentage = std_toeoff;
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).all_percentages = all_toeoff_percentages;
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).num_cycles = length(all_toeoff_percentages);
                    
                    fprintf('      Toe-off moyen pour %s: %.2f%% ± %.2f%% (%d cycles)\n', ...
                            Condition{iC}, mean_toeoff, std_toeoff, length(all_toeoff_percentages));
                end

                % Identification des outliers globaux pour Plat
                fprintf('      Identification des outliers globaux Plat...\n');
                
                all_outlier_indices = [];
                for m = 1:length(muscles)
                    muscle = muscles{m};
                    if isfield(CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}), muscle)
                        outlier_indices = CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle);
                        all_outlier_indices = union(all_outlier_indices, find(outlier_indices));
                    end
                end
                
                % Calcul des cycles moyens bruts et des maxima de référence pour Plat
                fprintf('      Calcul des cycles moyens bruts et maxima de référence...\n');
                
                for m = 1:length(muscles)
                    muscle = muscles{m};
                    
                    if isfield(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}), muscle)
                        all_cycles = CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle);
                        
                        % Déterminer les cycles valides
                        if ~isempty(all_outlier_indices)
                            total_cycles_available = length(all_cycles);
                            valid_cycle_indices = setdiff(1:total_cycles_available, all_outlier_indices);
                        else
                            valid_cycle_indices = 1:length(all_cycles);
                        end
                        
                        if ~isempty(valid_cycle_indices)
                            % Calcul du cycle moyen brut
                            all_points = zeros(length(valid_cycle_indices), 100);
                            for cycle_idx = 1:length(valid_cycle_indices)
                                cycle_num = valid_cycle_indices(cycle_idx);
                                if cycle_num <= length(all_cycles)
                                    all_points(cycle_idx, :) = all_cycles{cycle_num};
                                end
                            end
                            
                            cycle_moyen_brut = mean(all_points, 1);
                            max_cycle_moyen_brut = max(cycle_moyen_brut);
                            
                            % Stockage cycle moyen brut
                            if ~isfield(CYCLES_MOYENS_BRUTS, Participant{iP})
                                CYCLES_MOYENS_BRUTS.(Participant{iP}) = struct();
                            end
                            if ~isfield(CYCLES_MOYENS_BRUTS.(Participant{iP}), Condition{iC})
                                CYCLES_MOYENS_BRUTS.(Participant{iP}).(Condition{iC}) = struct();
                            end
                            if ~isfield(CYCLES_MOYENS_BRUTS.(Participant{iP}).(Condition{iC}), jambes{j})
                                CYCLES_MOYENS_BRUTS.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
                            end
                            
                            CYCLES_MOYENS_BRUTS.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = cycle_moyen_brut;
                            
                            % Stockage max de référence pour Plat
                            if ~isfield(MAX_CYCLES_BRUTS, Participant{iP})
                                MAX_CYCLES_BRUTS.(Participant{iP}) = struct();
                            end
                            if ~isfield(MAX_CYCLES_BRUTS.(Participant{iP}), jambes{j})
                                MAX_CYCLES_BRUTS.(Participant{iP}).(jambes{j}) = struct();
                            end
                            
                            MAX_CYCLES_BRUTS.(Participant{iP}).(jambes{j}).(muscle) = max_cycle_moyen_brut;
                            
                            fprintf('        Muscle %s: Max de référence = %.4f\n', muscle, max_cycle_moyen_brut);
                        end
                    end
                end
                
                break; % Sortir de la boucle après avoir traité Plat
            end
        end
        
        % ÉTAPE 2: Traiter toutes les conditions avec la normalisation appropriée
        fprintf('  === ÉTAPE 2: Traitement de toutes les conditions avec normalisation différentielle ===\n');
        
        for iC = 1:length(Condition)
            fprintf('    Condition: %s\n', Condition{iC});
            
            % Dossier d'enregistrement
            participant_folder = fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', Participant{iP}, Condition{iC});
            if ~exist(participant_folder, 'dir')
                mkdir(participant_folder);
            end

            % Si ce n'est pas Plat, extraire les cycles pour cette condition
            if ~strcmp(Condition{iC}, 'Plat')
                fprintf('      Extraction des cycles pour %s...\n', Condition{iC});
                fprintf('      Extraction des toe-offs %s...\n', Condition{iC});
                all_toeoff_percentages = [];
                
                % Traitement muscle par muscle
                for m = 1:length(muscles)
                    muscle = muscles{m};
                    
                    % Créer figure pour visualisation
                    figure_concatenated = figure('Name', sprintf('Cycles Concaténés - %s - %s - %s - %s', Participant{iP}, Condition{iC}, jambes{j}, muscle));
                    hold on;
                    
                    cycle_offset = 0;
                    cycle_positions = [];
                    muscle_maxima = [];
                    total_cycles = 0;
                    all_cycles_data = {}; % Cell array pour stocker tous les cycles

                    % Boucle des essais
                    for iEs = 1:length(Essai)  
                        file = [Participant{iP} '_' Condition{iC} '_' Essai{iEs} '.c3d'];

                        try
                            % EXTRACTION DES DONNEES BRUTES
                            data = btkReadAcquisition(file);
                            analogs = btkGetAnalogs(data);
                            sensor_name = sensor_association.(muscle);
                            EMG_signal = analogs.(sensor_name);

                            % Traitement du signal
                            FreqS = btkGetAnalogFrequency(data);
                            FreqVicon = 100;
                            EMGfilt = filtrage(EMG_signal, FreqS, 20, 450);
                            artifacts_info = characterizeArtifacts(EMGfilt,FreqS);
                            [signal_cleaned, EMGenvcleaned] = cleanEMGdata(EMGfilt, FreqS, 'STD', 5, 15, false);

                            % Détection HS
                            if strcmp(jambes{j}, 'left')
                                HS = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_signal);
                                TO = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                            else
                                HS = indiceRight(data, analogs, FreqS, FreqVicon, EMG_signal);
                                TO = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_signal); % TOE-OFF
                            end
                            
                            num_cycles = length(HS) - 1;
                            cycles_frames = cell(num_cycles, 1);

                            % === CALCUL DES POURCENTAGES TOE-OFF POUR CET ESSAI ===
                            if m == 1 % Ne calculer qu'une seule fois par essai (pour le premier muscle)
                                for i = 1:num_cycles
                                    cycle_start = HS(i);
                                    cycle_end = HS(i + 1);
                                    cycle_duration = cycle_end - cycle_start;
                                    
                                    % Trouver le TO dans ce cycle
                                    TO_in_cycle = TO(TO > cycle_start & TO < cycle_end);
                                    
                                    if ~isempty(TO_in_cycle)
                                        % Prendre le premier TO dans le cycle
                                        TO_frame = TO_in_cycle(1);
                                        % Calculer le pourcentage
                                        TO_percentage = ((TO_frame - cycle_start) / cycle_duration) * 100;
                                        all_toeoff_percentages = [all_toeoff_percentages, TO_percentage];
                                    else
                                        % Si pas de TO trouvé, utiliser une valeur par défaut
                                        warning('Pas de toe-off trouvé pour le cycle %d, essai %s, condition %s', i, Essai{iEs}, Condition{iC});
                                        all_toeoff_percentages = [all_toeoff_percentages, 60]; % Valeur typique
                                    end
                                end
                            end

                            for i = 1:num_cycles
                                cycles_frames{i} = HS(i):HS(i + 1);
                            end

                            % Traitement des cycles individuels
                            for i = 1:num_cycles
                                cycle_signal = EMGenvcleaned(cycles_frames{i});
                                cycle_signal_interpolated = interp1(linspace(1, length(cycle_signal), length(cycle_signal)), cycle_signal, linspace(1, length(cycle_signal), 100), 'pchip');
                                
                                % Stocker le cycle BRUT (non normalisé)
                                total_cycles = total_cycles + 1;
                                all_cycles_data{total_cycles} = cycle_signal_interpolated;
                                
                                cycle_positions = [cycle_positions, cycle_offset + length(cycle_signal_interpolated)];
                                muscle_maxima = [muscle_maxima, max(cycle_signal_interpolated)];

                                % Tracé individuel de chaque cycle
                                plot(cycle_offset + (1:length(cycle_signal_interpolated)), cycle_signal_interpolated, 'b');
                                cycle_offset = cycle_offset + length(cycle_signal_interpolated);
                            end
                            
                        catch ME
                            warning('Erreur lors du traitement du fichier %s: %s', file, ME.message);
                            continue;
                        end
                    end

                    % Finaliser la figure
                    y_limits = ylim;
                    for i = 1:length(cycle_positions)
                        plot([cycle_positions(i), cycle_positions(i)], y_limits, 'r', 'LineWidth', 0.5);
                    end

                    title_str = sprintf('All cycles for %s - %s - Muscle: %s - Leg: %s', ...
                                        Participant{iP}, Condition{iC}, muscle, jambes{j});
                    title(title_str, 'Interpreter', 'none');
                    xlabel('Temps normalisé');
                    ylabel('Enveloppe');

                    % Détection outliers basée sur les maxima
                    outlier_indices = detect_outliers(muscle_maxima, cycle_offset, 3);

                    % Sauvegarde figure
                    save_path = fullfile(participant_folder, sprintf('%s_%s_%s_%s.png', Participant{iP}, Condition{iC}, jambes{j}, muscle));
                    print(figure_concatenated, save_path, '-dpng', '-r150');
                    close(figure_concatenated);

                    % Stockage des données et outliers
                    if ~isfield(CYCLES_SIGNAL, Participant{iP})
                        CYCLES_SIGNAL.(Participant{iP}) = struct();
                    end
                    if ~isfield(CYCLES_SIGNAL.(Participant{iP}), Condition{iC})
                        CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}) = struct();
                    end
                    if ~isfield(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}), jambes{j})
                        CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
                    end
                    
                    CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = all_cycles_data;
                    
                    % Structures outliers
                    if ~isfield(CYCLES_OUTLIERS, Participant{iP})
                        CYCLES_OUTLIERS.(Participant{iP}) = struct();
                    end
                    if ~isfield(CYCLES_OUTLIERS.(Participant{iP}), Condition{iC})
                        CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}) = struct();
                    end
                    if ~isfield(CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}), jambes{j})
                        CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
                    end
                    
                    CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = outlier_indices;
                end

                % === STOCKAGE TOE-OFF MOYEN POUR CETTE CONDITION ===
                if ~isempty(all_toeoff_percentages)
                    mean_toeoff = mean(all_toeoff_percentages);
                    std_toeoff = std(all_toeoff_percentages);
                    
                    % Créer la structure CYCLES_TOEOFF
                    if ~isfield(CYCLES_TOEOFF, Participant{iP})
                        CYCLES_TOEOFF.(Participant{iP}) = struct();
                    end
                    if ~isfield(CYCLES_TOEOFF.(Participant{iP}), Condition{iC})
                        CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}) = struct();
                    end
                    
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).mean_percentage = mean_toeoff;
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).std_percentage = std_toeoff;
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).all_percentages = all_toeoff_percentages;
                    CYCLES_TOEOFF.(Participant{iP}).(Condition{iC}).(jambes{j}).num_cycles = length(all_toeoff_percentages);
                    
                    fprintf('      Toe-off moyen pour %s: %.2f%% ± %.2f%% (%d cycles)\n', ...
                            Condition{iC}, mean_toeoff, std_toeoff, length(all_toeoff_percentages));
                end
            end

            % Construction de la matrice de synergies avec normalisation différentielle
            fprintf('      Construction matrice synergies avec normalisation différentielle...\n');
            
            % Identification des outliers globaux
            all_outlier_indices = [];
            for m = 1:length(muscles)
                muscle = muscles{m};
                if isfield(CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}), muscle)
                    outlier_indices = CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle);
                    all_outlier_indices = union(all_outlier_indices, find(outlier_indices));
                end
            end
            
            % Déterminer le nombre de cycles valides
            if ~isempty(all_outlier_indices) && isfield(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}), muscles{1})
                total_cycles_available = length(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscles{1}));
                valid_cycle_indices = setdiff(1:total_cycles_available, all_outlier_indices);
            else
                if isfield(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}), muscles{1})
                    total_cycles_available = length(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscles{1}));
                    valid_cycle_indices = 1:total_cycles_available;
                else
                    valid_cycle_indices = [];
                end
            end
            
            num_valid_cycles = length(valid_cycle_indices);
            fprintf('        Nombre de cycles valides pour %s: %d\n', Condition{iC}, num_valid_cycles);
            
            if num_valid_cycles > 0
                % Préparer la matrice de synergies
                synergy_matrix = zeros(num_valid_cycles * 100, length(muscles));
                
                for m = 1:length(muscles)
                    muscle = muscles{m};
                    
                    % Récupérer tous les cycles de ce muscle
                    if isfield(CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}), muscle)
                        all_cycles = CYCLES_SIGNAL.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle);
                        
                        % Construire la colonne pour ce muscle avec normalisation différentielle
                        muscle_column = [];
                        
                        % NORMALISATION DIFFÉRENTIELLE SELON LA CONDITION
                        if strcmp(Condition{iC}, 'Plat')
                            % CONDITION PLAT: Normalisation par le maximum de chaque cycle individuel
                            fprintf('          Normalisation Plat: par maximum individuel de chaque cycle\n');
                            
                            for cycle_idx = 1:length(valid_cycle_indices)
                                cycle_num = valid_cycle_indices(cycle_idx);
                                if cycle_num <= length(all_cycles)
                                    cycle_data = all_cycles{cycle_num};
                                    
                                    % Normalisation par le maximum du cycle
                                    cycle_max = max(cycle_data);
                                    if cycle_max > 0
                                        cycle_normalized = cycle_data / cycle_max;
                                    else
                                        cycle_normalized = cycle_data;
                                        warning('Cycle avec maximum nul détecté pour %s - %s - %s - %s, cycle %d', ...
                                                Participant{iP}, Condition{iC}, jambes{j}, muscle, cycle_num);
                                    end
                                    
                                    % Concaténer ce cycle normalisé
                                    muscle_column = [muscle_column; cycle_normalized(:)];
                                end
                            end
                            
                        else
                            % CONDITIONS MEDIUM ET HIGH: Normalisation par le maximum du cycle moyen brut de Plat
                            fprintf('          Normalisation %s: par maximum cycle moyen Plat de référence\n', Condition{iC});
                            
                            % Récupérer le maximum de référence de la condition Plat
                            if isfield(MAX_CYCLES_BRUTS, Participant{iP}) && ...
                               isfield(MAX_CYCLES_BRUTS.(Participant{iP}), jambes{j}) && ...
                               isfield(MAX_CYCLES_BRUTS.(Participant{iP}).(jambes{j}), muscle)
                                
                                max_plat_reference = MAX_CYCLES_BRUTS.(Participant{iP}).(jambes{j}).(muscle);
                                fprintf('            Max référence Plat pour %s: %.4f\n', muscle, max_plat_reference);
                                
                                for cycle_idx = 1:length(valid_cycle_indices)
                                    cycle_num = valid_cycle_indices(cycle_idx);
                                    if cycle_num <= length(all_cycles)
                                        cycle_data = all_cycles{cycle_num}; % Cycle brut non normalisé
                                        
                                        % Normalisation par le maximum de référence Plat
                                        if max_plat_reference > 0
                                            cycle_normalized = cycle_data / max_plat_reference;
                                        else
                                            cycle_normalized = cycle_data;
                                            warning('Maximum de référence Plat nul pour %s - %s - %s', ...
                                                    Participant{iP}, jambes{j}, muscle);
                                        end
                                        
                                        % Concaténer ce cycle normalisé
                                        muscle_column = [muscle_column; cycle_normalized(:)];
                                    end
                                end
                                
                            else
                                warning('Pas de maximum de référence Plat trouvé pour %s - %s - %s. Utilisation normalisation individuelle par défaut.', ...
                                        Participant{iP}, jambes{j}, muscle);
                                
                                % Fallback: normalisation individuelle
                                for cycle_idx = 1:length(valid_cycle_indices)
                                    cycle_num = valid_cycle_indices(cycle_idx);
                                    if cycle_num <= length(all_cycles)
                                        cycle_data = all_cycles{cycle_num};
                                        cycle_max = max(cycle_data);
                                        if cycle_max > 0
                                            cycle_normalized = cycle_data / cycle_max;
                                        else
                                            cycle_normalized = cycle_data;
                                        end
                                        muscle_column = [muscle_column; cycle_normalized(:)];
                                    end
                                end
                            end
                        end
                        
                        % Stocker dans la matrice de synergies
                        if length(muscle_column) == num_valid_cycles * 100
                            synergy_matrix(:, m) = muscle_column;
                        else
                            warning('Taille incorrecte pour le muscle %s: attendu %d, obtenu %d', ...
                                    muscle, num_valid_cycles * 100, length(muscle_column));
                        end
                        
                        % === Calcul du cycle moyen à partir de la matrice de synergies ===
n_cycles = num_valid_cycles;
n_points = 100;

muscle_vector = synergy_matrix(:, m); % Extraire les données concaténées du muscle m

if length(muscle_vector) == n_cycles * n_points
    % Reshape en [n_cycles x 100]
    muscle_matrix = reshape(muscle_vector, [n_points, n_cycles])';

    % Calcul du cycle moyen
    mean_cycle = mean(muscle_matrix, 1);

    % Stockage dans la structure CYCLES_MOYENS
    if ~isfield(CYCLES_MOYENS, Participant{iP})
        CYCLES_MOYENS.(Participant{iP}) = struct();
    end
    if ~isfield(CYCLES_MOYENS.(Participant{iP}), Condition{iC})
        CYCLES_MOYENS.(Participant{iP}).(Condition{iC}) = struct();
    end
    if ~isfield(CYCLES_MOYENS.(Participant{iP}).(Condition{iC}), jambes{j})
        CYCLES_MOYENS.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
    end

    CYCLES_MOYENS.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = mean_cycle;

% === Sauvegarde figure du cycle moyen + SD ===
    x = linspace(0, 100, n_points);
    sd_cycle = std(muscle_matrix, 0, 1);

    fig = figure('Visible', 'off'); % ne pas afficher
    hold on;

    % Zone écart-type
    fill([x, fliplr(x)], ...
         [mean_cycle + sd_cycle, fliplr(mean_cycle - sd_cycle)], ...
         [1 0 0], 'FaceAlpha', 0.3, 'EdgeColor', 'none');

    % Courbe moyenne
    plot(x, mean_cycle, 'r', 'LineWidth', 2);

    title(sprintf('Cycle moyen normalisé - %s - %s - Muscle: %s - Leg: %s', ...
          Participant{iP}, Condition{iC}, muscle, jambes{j}), 'Interpreter', 'none');
    xlabel('% du cycle');
    ylabel('Activation normalisée');
    xlim([0 100]);

    % Créer dossier s’il n’existe pas
    fig_folder = fullfile(participant_folder, 'CycleMoyen');
    if ~exist(fig_folder, 'dir')
        mkdir(fig_folder);
    end

    % Sauvegarde
    fig_path = fullfile(fig_folder, sprintf('%s_%s_%s_%s_cycleMoyen.png', ...
                    Participant{iP}, Condition{iC}, jambes{j}, muscle));
    print(fig, fig_path, '-dpng', '-r150');
    close(fig);

else
    warning('La taille du vecteur de données ne correspond pas au nombre attendu de points pour %s - %s - %s - %s', ...
            Participant{iP}, Condition{iC}, jambes{j}, muscle);
end

                        % Stockage des cycles normalisés avec la nouvelle méthode
                        if ~isfield(CYCLES_SIGNAL_NORMALIZED, Participant{iP})
                            CYCLES_SIGNAL_NORMALIZED.(Participant{iP}) = struct();
                        end
                        if ~isfield(CYCLES_SIGNAL_NORMALIZED.(Participant{iP}), Condition{iC})
                            CYCLES_SIGNAL_NORMALIZED.(Participant{iP}).(Condition{iC}) = struct();
                        end
                        if ~isfield(CYCLES_SIGNAL_NORMALIZED.(Participant{iP}).(Condition{iC}), jambes{j})
                            CYCLES_SIGNAL_NORMALIZED.(Participant{iP}).(Condition{iC}).(jambes{j}) = struct();
                        end
                        
                     % Stocker les cycles normalisés selon la nouvelle méthode
                        normalized_cycles = {};
                        if strcmp(Condition{iC}, 'Plat')
                            % Plat: normalisation individuelle par le max de chaque cycle
                            for cycle_idx = 1:length(valid_cycle_indices)
                                cycle_num = valid_cycle_indices(cycle_idx);
                                if cycle_num <= length(all_cycles)
                                    cycle_data = all_cycles{cycle_num};
                                    cycle_max = max(cycle_data);
                                    if cycle_max > 0
                                        normalized_cycles{cycle_idx} = cycle_data / cycle_max;
                                    else
                                        normalized_cycles{cycle_idx} = cycle_data;
                                    end
                                end
                            end
                        else
                            % Medium et High: normalisation POINT PAR POINT par le cycle moyen Plat
                            if isfield(CYCLES_MOYENS_BRUTS, Participant{iP}) && ...
                               isfield(CYCLES_MOYENS_BRUTS.(Participant{iP}), 'Plat') && ...
                               isfield(CYCLES_MOYENS_BRUTS.(Participant{iP}).('Plat'), jambes{j}) && ...
                               isfield(CYCLES_MOYENS_BRUTS.(Participant{iP}).('Plat').(jambes{j}), muscle)
                                
                                cycle_plat_reference = CYCLES_MOYENS_BRUTS.(Participant{iP}).('Plat').(jambes{j}).(muscle);
                                
                                for cycle_idx = 1:length(valid_cycle_indices)
                                    cycle_num = valid_cycle_indices(cycle_idx);
                                    if cycle_num <= length(all_cycles)
                                        cycle_data = all_cycles{cycle_num}; % Cycle brut
                                        
                                        % Normalisation POINT PAR POINT
                                        if length(cycle_data) == length(cycle_plat_reference)
                                            cycle_normalized = zeros(size(cycle_data));
                                            for pt = 1:length(cycle_data)
                                                if cycle_plat_reference(pt) > 0
                                                    cycle_normalized(pt) = cycle_data(pt) / cycle_plat_reference(pt);
                                                else
                                                    cycle_normalized(pt) = cycle_data(pt);
                                                end
                                            end
                                            normalized_cycles{cycle_idx} = cycle_normalized;
                                        else
                                            normalized_cycles{cycle_idx} = cycle_data;
                                        end
                                    end
                                end
                            end
                        end
                        
                        CYCLES_SIGNAL_NORMALIZED.(Participant{iP}).(Condition{iC}).(jambes{j}).(muscle) = normalized_cycles;
                      
                    else
                        warning('Pas de données de cycles trouvées pour %s - %s - %s - %s', ...
                                Participant{iP}, Condition{iC}, jambes{j}, muscle);
                    end
                end
                
                % Stockage de la matrice de synergies
                if ~isfield(SYNERGY_MATRIX, Participant{iP})
                    SYNERGY_MATRIX.(Participant{iP}) = struct();
                end
                if ~isfield(SYNERGY_MATRIX.(Participant{iP}), Condition{iC})
                    SYNERGY_MATRIX.(Participant{iP}).(Condition{iC}) = struct();
                end
                
                SYNERGY_MATRIX.(Participant{iP}).(Condition{iC}).(jambes{j}) = synergy_matrix;
                
                % Stockage du nombre de cycles
                if ~isfield(CYCLES_COUNT, Participant{iP})
                    CYCLES_COUNT.(Participant{iP}) = struct();
                end
                if ~isfield(CYCLES_COUNT.(Participant{iP}), Condition{iC})
                    CYCLES_COUNT.(Participant{iP}).(Condition{iC}) = struct();
                end
                
                CYCLES_COUNT.(Participant{iP}).(Condition{iC}).(jambes{j}) = num_valid_cycles;
                
                fprintf('        Matrice de synergies créée: %d x %d (cycles x muscles)\n', ...
                        size(synergy_matrix, 1), size(synergy_matrix, 2));
                
            else
                warning('Aucun cycle valide trouvé pour %s - %s - %s', ...
                        Participant{iP}, Condition{iC}, jambes{j});
            end
        end
    end
end

% Restaurer l'affichage des figures
set(0, 'DefaultFigureVisible', original_visible);

fprintf('\n=== SAUVEGARDE DES RÉSULTATS ===\n');

% Définir le chemin de sauvegarde
save_path = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS';
if ~exist(save_path, 'dir')
    mkdir(save_path);
end

% Sauvegarde des structures principales
save_file = fullfile(save_path, [Participant{iP} '_MATRIX.mat']);
save(save_file, 'SYNERGY_MATRIX', 'CYCLES_SIGNAL_NORMALIZED', ...
     'CYCLES_MOYENS_BRUTS', 'CYCLES_MOYENS', 'CYCLES_COUNT', 'CYCLES_OUTLIERS', 'CYCLES_TOEOFF', '-v7.3');

fprintf('Structures sauvegardées dans: %s\n', save_file);

fprintf('\n=== TRAITEMENT TERMINÉ ===\n');
fprintf('Script terminé avec succès!\n');
fprintf('Normalisation différentielle appliquée:\n');
fprintf('- Plat: normalisation par max individuel de chaque cycle\n');
fprintf('- Medium/High: normalisation POINT PAR POINT par cycle moyen brut Plat de référence\n');