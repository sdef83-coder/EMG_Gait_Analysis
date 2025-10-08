%% DRESSER LES CYCLES MOYENS DES PARTICIPANTS EN DISSOCIANT LES JAMBES
% Affiche le nom du participant avec la valeur max à chaque 10%
% Attention changer ligne 9, 148, 186

clc; clear; close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix')

% Liste des participants, conditions, muscles et jambes
Participant = {'CTL_14', 'CTL_15', 'CTL_23', 'CTL_37', 'CTL_40', 'CTL_53', 'CTL_63'}; 
% Adultes : 'CTL_04', 'CTL_05', 'CTL_07', 'CTL_08', 'CTL_12', 'CTL_18', 'CTL_24', 'CTL_25', 'CTL_27', 'CTL_28', 'CTL_32'
% Enfants : 'CTL_01', 'CTL_02', 'CTL_06', 'CTL_09', 'CTL_10', 'CTL_11', 'CTL_16', 'CTL_19', 'CTL_20', 'CTL_34'
% Jeunes Enfants : 'CTL_14', 'CTL_15', 'CTL_23', 'CTL_37', 'CTL_40'
Conditions = {'Plat', 'Medium', 'High'};  
Muscles = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'};
Jambes = {'left', 'right'};  

% Initialisation des structures
mean_cycles_global = struct();
participant_tracking = struct();

% Charger les données de tous les participants
for iP = 1:length(Participant)
    filename = [Participant{iP} '_MATRIX', '.mat'];
    metadonnees = load(filename);
    
    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        for j = 1:length(Jambes)
            jambe = Jambes{j};
            for iM = 1:length(Muscles)
                muscle = Muscles{iM};

                if isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe), muscle)
                    cycle_mean = metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe).(muscle);

                    % Initialiser les structures si nécessaire
                    if ~isfield(mean_cycles_global, condition)
                        mean_cycles_global.(condition) = struct();
                        participant_tracking.(condition) = struct();
                    end
                    if ~isfield(mean_cycles_global.(condition), muscle)
                        mean_cycles_global.(condition).(muscle) = struct('left', [], 'right', []);
                        participant_tracking.(condition).(muscle) = struct('left', {{}}, 'right', {{}});
                    end
                    if ~isfield(mean_cycles_global.(condition).(muscle), jambe)
                        mean_cycles_global.(condition).(muscle).(jambe) = [];
                        participant_tracking.(condition).(muscle).(jambe) = {};
                    end

                    % Stocker les données
                    mean_cycles_global.(condition).(muscle).(jambe) = [mean_cycles_global.(condition).(muscle).(jambe); cycle_mean];
                    participant_tracking.(condition).(muscle).(jambe){end+1} = Participant{iP};
                else
                    fprintf('Le muscle %s est absent pour le participant %s dans la condition %s et la jambe %s.\n', muscle, Participant{iP}, condition, jambe);
                end
            end
        end
    end
end

% Visualisation avec labels des participants ayant les valeurs max à chaque 10%
for iM = 1:length(Muscles)
    muscle = Muscles{iM};

    for j = 1:length(Jambes)
        jambe = Jambes{j};

        figure;
        hold on;

        colors = {[0 0 1], [0 0.5 0], [1 0 0]};  % bleu, vert, rouge
        legends = {};
        h_legend = [];

        for iC = 1:length(Conditions)
            condition = Conditions{iC};

            if isfield(mean_cycles_global.(condition).(muscle), jambe)
                all_participant_cycles = mean_cycles_global.(condition).(muscle).(jambe);
                participant_names = participant_tracking.(condition).(muscle).(jambe);
                
                valid_rows = all(isfinite(all_participant_cycles), 2);
                all_participant_cycles = all_participant_cycles(valid_rows, :);
                participant_names = participant_names(valid_rows);

                % Afficher le nombre de cycles valides/exclus
                total_cycles = size(mean_cycles_global.(condition).(muscle).(jambe), 1);
                valid_cycles = size(all_participant_cycles, 1);
                fprintf('[%s - %s - %s] %d cycles valides sur %d (exclus : %d)\n', ...
                    condition, muscle, jambe, valid_cycles, total_cycles, total_cycles - valid_cycles);

                mean_cycle_all_participants = mean(all_participant_cycles, 1, 'omitnan');
                std_cycle_condition = std(all_participant_cycles, 0, 1, 'omitnan');

                % Tracer les cycles individuels (courbes fines, transparentes)
                for k = 1:size(all_participant_cycles, 1)
                    plot(all_participant_cycles(k, :), 'Color', [colors{iC} 0.3], 'LineWidth', 0.5);
                end

                % Identifier et afficher les participants avec valeurs max à chaque 10%
                cycle_length = size(all_participant_cycles, 2);
                step_size = round(cycle_length / 10);  % Diviser en 10 segments
                
                for step = 1:10
                    % Position dans le cycle (en %)
                    position_index = step * step_size;
                    if position_index > cycle_length
                        position_index = cycle_length;
                    end
                    
                    % Trouver le participant avec la valeur maximale à cette position
                    [max_value, max_participant_idx] = max(all_participant_cycles(:, position_index));
                    max_participant_name = participant_names{max_participant_idx};
                    
                    % Afficher le nom du participant à cette position
                    text(position_index, max_value + 0.05, max_participant_name, ...
                         'FontSize', 8, 'Color', colors{iC}, 'FontWeight', 'bold', ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                         'BackgroundColor', 'white', 'EdgeColor', colors{iC}, ...
                         'Margin', 1);
                    
                    % Marquer le point avec un cercle
                    scatter(position_index, max_value, 40, colors{iC}, 'filled', ...
                           'MarkerEdgeColor', 'white', 'LineWidth', 1);   
                end

                % Tracer la courbe de la moyenne pour cette condition et cette jambe
                h = plot(mean_cycle_all_participants, 'LineWidth', 2, 'Color', colors{iC}, 'DisplayName', condition);

                % Afficher la zone d'écart type autour de la moyenne
                lower_bound = max(0, mean_cycle_all_participants - std_cycle_condition);
                fill([1:length(mean_cycle_all_participants), fliplr(1:length(mean_cycle_all_participants))], ...
                    [mean_cycle_all_participants + std_cycle_condition, fliplr(lower_bound)], ...
                    colors{iC}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');

                % Compter les participants uniques (pas les cycles)
                unique_participants = unique(participant_names);
                nb_participants_uniques = length(unique_participants);
                
                % Ajouter à la légende avec le nombre de participants uniques
                legends{end+1} = [condition ' (n=' num2str(nb_participants_uniques) ')'];
                h_legend(end+1) = h;
            end
        end

        % Ajouter des légendes et un titre
        legend(h_legend, legends, 'Location', 'best');
        title(['Cycle moyen ' muscle ' - ' jambe ' - Population Adulte']);
        xlabel('Temps normalisé');
        ylabel('EMG normalisé');
        grid on;

        % Calculer le nombre total de participants uniques pour ce muscle et cette jambe
        all_unique_participants = {};
        for iC = 1:length(Conditions)
            condition = Conditions{iC};
            if isfield(mean_cycles_global.(condition).(muscle), jambe)
                all_participant_cycles = mean_cycles_global.(condition).(muscle).(jambe);
                participant_names = participant_tracking.(condition).(muscle).(jambe);
                valid_rows = all(isfinite(all_participant_cycles), 2);
                valid_participant_names = participant_names(valid_rows);
                all_unique_participants = [all_unique_participants, valid_participant_names];
            end
        end
        total_unique_participants = length(unique(all_unique_participants));
        
        % Afficher le nombre total de participants uniques en bas à droite
        ax = gca;
        xlim_val = xlim;
        ylim_val = ylim;
        text(xlim_val(2)*0.95, ylim_val(1)*0.95 + ylim_val(2)*0.05, ...
             ['Participants uniques: ' num2str(total_unique_participants)], ...
             'FontSize', 9, 'FontWeight', 'bold', ...
             'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', ...
             'BackgroundColor', 'white', 'EdgeColor', 'black', ...
             'Margin', 3);

        hold off;

        % Enregistrer la figure dans le dossier de résultats
        saveas(gcf, fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ['Cycle_Moyen_' muscle '_' jambe '_MaxLabels.png']));
    end
end

% Sauvegarder les résultats finaux
save('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\mean_cycles_adultes.mat', 'mean_cycles_global');

%% COMBINE LES JAMBES - APERCU DE L'ACTIVATION

clc; clear; close all;
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\FILTERED');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script'));

% === CONFIGURATION DU GROUPE À ÉTUDIER ===
groupe_a_etudier = 'JeunesEnfants';  % Modifier ici : 'Adultes', 'Adolescents', 'Enfants', 'JeunesEnfants'

% Charger la structure Group définie dans ParticipantGroup.m
ParticipantGroup;

% Extraire la liste des participants correspondant au groupe sélectionné
Participant = Group.(groupe_a_etudier);

% Affichage
fprintf('\n=== Groupe sélectionné : %s ===\n', groupe_a_etudier);
fprintf('Participants inclus : %s\n\n', strjoin(Participant, ', '));

Conditions = {'Plat', 'Medium', 'High'};  
Muscles = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'};
Jambes = {'left', 'right'};  

% Initialisation d'une matrice pour stocker les cycles moyens par condition et muscle (sans distinction de jambe)
mean_cycles_combined = struct();
% Structure pour stocker les données par participant avec suivi des noms
participant_tracking = struct();

% Charger les données de tous les participants et calculer les cycles moyens par condition et muscle
for iP = 1:length(Participant)
    filename = [Participant{iP} '_MATRIX', '.mat'];  % Nom du fichier pour chaque participant
    metadonnees = load(filename);  % Charger le fichier MATLAB
    
    % Parcourir chaque condition
    for iC = 1:length(Conditions)
        condition = Conditions{iC};

        % Parcourir chaque muscle
        for iM = 1:length(Muscles)
            muscle = Muscles{iM};
            
            % Initialiser un tableau pour stocker les cycles des deux jambes pour ce participant
            participant_cycles = [];
            valid_cycles = []; % Pour stocker uniquement les cycles valides (sans NaN)

            % Parcourir chaque jambe pour collecter les données
            for j = 1:length(Jambes)
                jambe = Jambes{j};

                % Vérifier si le champ existe pour ce participant, condition, jambe et muscle
                if isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe), muscle)
                    % Récupérer le cycle moyen normalisé pour ce muscle, cette condition et cette jambe
                    cycle_mean = metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe).(muscle);
                    
                    % Vérifier si le cycle contient des NaN ou INF
                   if all(isfinite(cycle_mean))
                        % Si pas de NaN, ajouter à la liste des cycles valides
                        valid_cycles = [valid_cycles; cycle_mean];
                        fprintf('Cycle valide ajouté pour %s - %s - %s - %s\n', Participant{iP}, condition, jambe, muscle);
                    else
                        fprintf('Cycle avec NaN ignoré pour %s - %s - %s - %s\n', Participant{iP}, condition, jambe, muscle);
                    end
                else
                    fprintf('Le muscle %s est absent pour le participant %s dans la condition %s et la jambe %s.\n', muscle, Participant{iP}, condition, jambe);
                end
            end
            
            % Si au moins un cycle valide a été trouvé pour ce participant
            if ~isempty(valid_cycles)
                % Calculer la moyenne des cycles valides seulement
                participant_mean_cycle = mean(valid_cycles, 1);
                
                % Vérifier une dernière fois que le résultat ne contient pas de NaN
                if all(isfinite(participant_mean_cycle))
                    % Stocker les cycles moyens dans la structure globale
                    if ~isfield(mean_cycles_combined, condition)
                        mean_cycles_combined.(condition) = struct();  % Si la condition n'existe pas encore
                        participant_tracking.(condition) = struct();
                    end
                    if ~isfield(mean_cycles_combined.(condition), muscle)
                        mean_cycles_combined.(condition).(muscle) = [];  % Si le muscle n'existe pas encore
                        participant_tracking.(condition).(muscle) = {};
                    end
                    
                    % Ajouter le cycle moyen combiné à la liste pour ce muscle et condition
                    mean_cycles_combined.(condition).(muscle) = [mean_cycles_combined.(condition).(muscle); participant_mean_cycle];
                    % Ajouter le nom du participant pour le suivi
                    participant_tracking.(condition).(muscle){end+1} = Participant{iP};
                    
                    fprintf('Cycle moyen calculé et ajouté pour %s - %s - %s\n', Participant{iP}, condition, muscle);
                else
                    fprintf('Cycle moyen avec NaN rejeté pour %s - %s - %s\n', Participant{iP}, condition, muscle);
                end
            else
                fprintf('Aucun cycle valide trouvé pour %s - %s - %s (les deux jambes ont des NaN)\n', Participant{iP}, condition, muscle);
            end
        end
    end
end

% Calculer les cycles moyens globaux pour chaque muscle et condition
global_mean_cycles = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    global_mean_cycles.(condition) = struct();
    
    for iM = 1:length(Muscles)
        muscle = Muscles{iM};
        if isfield(mean_cycles_combined, condition) && isfield(mean_cycles_combined.(condition), muscle)
            % Vérifier qu'il y a des données pour ce muscle et condition
            if ~isempty(mean_cycles_combined.(condition).(muscle))
                % Moyenne des cycles moyens pour ce muscle et cette condition (toutes jambes confondues)
                temp_mean = mean(mean_cycles_combined.(condition).(muscle), 1);
                
                % Vérifier que la moyenne globale ne contient pas de NaN
                if ~any(isnan(temp_mean))
                    global_mean_cycles.(condition).(muscle) = temp_mean;
                    fprintf('Cycle moyen global calculé pour %s - %s\n', condition, muscle);
                else
                    fprintf('Cycle moyen global avec NaN rejeté pour %s - %s\n', condition, muscle);
                end
            else
                fprintf('Aucune donnée disponible pour %s - %s\n', condition, muscle);
            end
        end
    end
end

% Normaliser les cycles moyens par rapport à la condition Plat
normalized_mean_cycles = struct();

for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    
    % Vérifier si la condition Plat existe pour ce muscle
    if isfield(global_mean_cycles, 'Plat') && isfield(global_mean_cycles.Plat, muscle)
        % Trouver la valeur maximale dans la condition Plat pour ce muscle
        max_flat = max(global_mean_cycles.Plat.(muscle));
        
        % Vérifier que max_flat n'est pas NaN ou zéro
        if ~isnan(max_flat) && max_flat > 0
            % Normaliser chaque condition par rapport à cette valeur maximale
            for iC = 1:length(Conditions)
                condition = Conditions{iC};
                
                if isfield(global_mean_cycles, condition) && isfield(global_mean_cycles.(condition), muscle)
                    if ~isfield(normalized_mean_cycles, condition)
                        normalized_mean_cycles.(condition) = struct();
                    end
                    
                    % Normaliser le cycle moyen de cette condition par rapport au max de la condition Plat
                    normalized_mean_cycles.(condition).(muscle) = global_mean_cycles.(condition).(muscle) / max_flat;
                    fprintf('Normalisation réussie pour %s - %s\n', condition, muscle);
                end
            end
        else
            fprintf('Impossible de normaliser %s : max_flat invalide (%f)\n', muscle, max_flat);
        end
    else
        fprintf('Condition Plat manquante pour le muscle %s - normalisation impossible\n', muscle);
    end
end

% Normaliser les données par participant
normalized_participant_cycles = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    if isfield(mean_cycles_combined, condition)
        normalized_participant_cycles.(condition) = struct();
        
        for iM = 1:length(Muscles)
            muscle = Muscles{iM};
            if isfield(mean_cycles_combined.(condition), muscle)
                % Récupérer la valeur max_flat pour ce muscle
                if isfield(global_mean_cycles, 'Plat') && isfield(global_mean_cycles.Plat, muscle)
                    max_flat = max(global_mean_cycles.Plat.(muscle));
                    
                    if ~isnan(max_flat) && max_flat > 0
                        % Normaliser tous les cycles des participants pour cette condition/muscle
                        normalized_participant_cycles.(condition).(muscle) = mean_cycles_combined.(condition).(muscle) / max_flat;
                    end
                end
            end
        end
    end
end

% === CALCUL DU TOEOFF MOYEN PAR CONDITION ===
toeoff_moyen = struct();  % contiendra toe-off par condition

for iC = 1:length(Conditions)
    condition = Conditions{iC};
    toeoff_condition = [];

    for iP = 1:length(Participant)
        participant = Participant{iP};
        filename = [participant '_MATRIX.mat'];
        data = load(filename);

        try
            left_toeoff = data.CYCLES_TOEOFF.(participant).(condition).left.mean_percentage;
            right_toeoff = data.CYCLES_TOEOFF.(participant).(condition).right.mean_percentage;

            if isfinite(left_toeoff) && isfinite(right_toeoff)
                mean_participant_toeoff = mean([left_toeoff, right_toeoff]);
                toeoff_condition(end+1) = mean_participant_toeoff;
            end
        catch
            fprintf('Toe-Off non dispo pour %s - %s\n', participant, condition);
        end
    end

    if ~isempty(toeoff_condition)
        toeoff_moyen.(condition) = mean(toeoff_condition);
        fprintf('Toe-Off moyen %s = %.2f%%\n', condition, toeoff_moyen.(condition));
    else
        toeoff_moyen.(condition) = NaN;
    end
end

% NOUVELLE PARTIE : Visualiser les cycles par participant avec la méthode du deuxième code
fprintf('\n=== CRÉATION DES GRAPHIQUES PAR PARTICIPANT ===\n');
for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    
    % Vérifier s'il y a des données à afficher pour ce muscle
    has_data = false;
    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        if isfield(normalized_participant_cycles, condition) && isfield(normalized_participant_cycles.(condition), muscle)
            has_data = true;
            break;
        end
    end
    
    if ~has_data
        fprintf('Aucune donnée à afficher pour le muscle %s - graphique ignoré\n', muscle);
        continue;
    end
    
    % Créer une figure pour ce muscle
    figure;
    hold on;
    
    % Couleurs pour chaque condition
    colors = {[0 0 1], [0 0.5 0], [1 0 0]};  % Bleu pour Plat, Vert pour Medium, Rouge pour High
    legends = {};
    h_legend = [];
    
    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        
        if isfield(normalized_participant_cycles, condition) && isfield(normalized_participant_cycles.(condition), muscle)
            all_participant_cycles = normalized_participant_cycles.(condition).(muscle);
            participant_names = participant_tracking.(condition).(muscle);
            
            % Filtrer les cycles valides
            valid_rows = all(isfinite(all_participant_cycles), 2);
            all_participant_cycles = all_participant_cycles(valid_rows, :);
            participant_names = participant_names(valid_rows);
            
            % Afficher le nombre de cycles valides/exclus
            total_cycles = size(normalized_participant_cycles.(condition).(muscle), 1);
            valid_cycles = size(all_participant_cycles, 1);
            fprintf('[%s - %s] %d cycles valides sur %d (exclus : %d)\n', ...
                condition, muscle, valid_cycles, total_cycles, total_cycles - valid_cycles);
            
            if ~isempty(all_participant_cycles)
                mean_cycle_all_participants = mean(all_participant_cycles, 1, 'omitnan');
                std_cycle_condition = std(all_participant_cycles, 0, 1, 'omitnan');
                
                % Tracer les cycles individuels (courbes fines, transparentes)
                for k = 1:size(all_participant_cycles, 1)
                    plot(all_participant_cycles(k, :), 'Color', [colors{iC} 0.3], 'LineWidth', 0.5);
                end
                
                % Identifier et afficher les participants avec valeurs max à chaque 10%
                cycle_length = size(all_participant_cycles, 2);
                step_size = round(cycle_length / 10);  % Diviser en 10 segments
                
                for step = 1:10
                    % Position dans le cycle (en %)
                    position_index = step * step_size;
                    if position_index > cycle_length
                        position_index = cycle_length;
                    end
                    
                    % Trouver le participant avec la valeur maximale à cette position
                    [max_value, max_participant_idx] = max(all_participant_cycles(:, position_index));
                    max_participant_name = participant_names{max_participant_idx};
                    
                    % Afficher le nom du participant à cette position
                    text(position_index, max_value + 0.05, max_participant_name, ...
                         'FontSize', 8, 'Color', colors{iC}, 'FontWeight', 'bold', ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                         'BackgroundColor', 'white', 'EdgeColor', colors{iC}, ...
                         'Margin', 1);
                    
                    % Marquer le point avec un cercle
                    scatter(position_index, max_value, 40, colors{iC}, 'filled', ...
                           'MarkerEdgeColor', 'white', 'LineWidth', 1);   
                end
                
                % Tracer la courbe de la moyenne pour cette condition
                h = plot(mean_cycle_all_participants, 'LineWidth', 2, 'Color', colors{iC}, 'DisplayName', condition);
                
                % Afficher la zone d'écart type autour de la moyenne
                lower_bound = max(0, mean_cycle_all_participants - std_cycle_condition);
                fill([1:length(mean_cycle_all_participants), fliplr(1:length(mean_cycle_all_participants))], ...
                    [mean_cycle_all_participants + std_cycle_condition, fliplr(lower_bound)], ...
                    colors{iC}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
                
                % Ajouter à la légende
                legends{end+1} = condition;
                h_legend(end+1) = h;
            end
        end
    end
    
    % Ajouter des légendes et un titre
    if ~isempty(h_legend)
        legend(h_legend, legends, 'Location', 'best');
        title(['Cycle moyen ' muscle ' - Toutes jambes - Population ' groupe_a_etudier]);
        xlabel('Cycle de marche (%)');
        ylabel('Activation musculaire normalisée');
        grid on;
    end
    
    hold off;
    
    % Enregistrer la figure
    if ~isempty(h_legend)
        saveas(gcf, fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ['Participants_' muscle '_MaxLabels.png']));
        fprintf('Graphique par participant sauvegardé pour le muscle %s\n', muscle);
    else
        close(gcf);
        fprintf('Aucune courbe à afficher pour le muscle %s - figure fermée\n', muscle);
    end
end

% === Dossier de sauvegarde des figures selon le groupe étudié ===
base_output_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle';
output_dir = fullfile(base_output_dir, groupe_a_etudier);
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Visualiser les cycles normalisés pour chaque muscle (graphique global existant)
fprintf('\n=== CRÉATION DES GRAPHIQUES GLOBAUX ===\n');
colors = {'b', 'g', 'r'};  % Plat, Medium, High

for iC = 1:length(Conditions)
    condition = Conditions{iC};
    figure('Name', ['Cycle - ' condition], 'Position', [100 100 1800 900]);
    
    for iM = 1:length(Muscles)
        muscle = Muscles{iM};
        subplot(2, 4, iM); hold on;
        title(strrep(muscle, 'EMG_', ''), 'Interpreter', 'none');
        xlabel('Cycle (%)'); ylabel('Activation normalisée');
        ylim([0 2]);
        
        % LIGNE VERTICALE TOEOFF MOYEN
        if isfield(toeoff_moyen, condition) && ~isnan(toeoff_moyen.(condition))
            x_toeoff = toeoff_moyen.(condition);
            y_limits = ylim;
            plot([x_toeoff, x_toeoff], y_limits, '--', 'Color', colors{iC}, 'LineWidth', 1.5);
        end


        if isfield(normalized_mean_cycles.(condition), muscle)
            cycle = normalized_mean_cycles.(condition).(muscle);
            if isfield(global_mean_cycles.Plat, muscle)
                max_flat = max(global_mean_cycles.Plat.(muscle));
                if ~isnan(max_flat) && max_flat > 0
                    std_cycle = nanstd(mean_cycles_combined.(condition).(muscle), 0, 1) / max_flat;
                else
                    std_cycle = zeros(size(cycle));
                end
            else
                std_cycle = zeros(size(cycle));
            end
            plot(cycle, 'LineWidth', 2, 'Color', colors{iC});
            fill([1:length(cycle), fliplr(1:length(cycle))], ...
                 [cycle + std_cycle, fliplr(cycle - std_cycle)], ...
                 colors{iC}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end

        grid on;
    end

    saveas(gcf, fullfile(output_dir, ['Cycle_Normalisé_' condition '_AllMuscles.png']));
    fprintf('Figure sauvegardée : Cycle_Normalisé_%s_AllMuscles.png\n', condition);
end

% === FIGURE COMBINÉE AVEC LES 3 CONDITIONS ===
figure('Name', 'Cycle - 3 Conditions', 'Position', [100 100 1800 900]);
for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    subplot(2, 4, iM); hold on;
    title(strrep(muscle, 'EMG_', ''), 'Interpreter', 'none');
    xlabel('Cycle (%)'); ylabel('Activation normalisée');
    ylim([0 2]);

    h_legend = [];
    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        if isfield(normalized_mean_cycles, condition) && isfield(normalized_mean_cycles.(condition), muscle)
            cycle = normalized_mean_cycles.(condition).(muscle);
            if isfield(global_mean_cycles.Plat, muscle)
                max_flat = max(global_mean_cycles.Plat.(muscle));
                if ~isnan(max_flat) && max_flat > 0
                    std_cycle = nanstd(mean_cycles_combined.(condition).(muscle), 0, 1) / max_flat;
                else
                    std_cycle = zeros(size(cycle));
                end
            else
                std_cycle = zeros(size(cycle));
            end
            h = plot(cycle, 'LineWidth', 2, 'Color', colors{iC});
            fill([1:length(cycle), fliplr(1:length(cycle))], ...
                 [cycle + std_cycle, fliplr(cycle - std_cycle)], ...
                 colors{iC}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            h_legend(end+1) = h;

            % LIGNE VERTICALE TOEOFF MOYEN
if isfield(toeoff_moyen, condition) && ~isnan(toeoff_moyen.(condition))
    x_toeoff = toeoff_moyen.(condition);
    y_limits = ylim;
    plot([x_toeoff, x_toeoff], y_limits, '--', 'Color', colors{iC}, 'LineWidth', 1.5);
end

        end
    end
    grid on;
    if iM == 1
        legend(h_legend, Conditions, 'Location', 'best');
    end
end

saveas(gcf, fullfile(output_dir, 'Cycle_Normalisé_AllConditions_AllMuscles.png'));
fprintf('Figure globale sauvegardée : Cycle_Normalisé_AllConditions_AllMuscles.png\n');

% Pour afficher qu'un seul muscle
% for iM = 1:length(Muscles)
%     muscle = Muscles{iM};
% 
%     % Vérifier s'il y a des données à afficher pour ce muscle
%     has_data = false;
%     for iC = 1:length(Conditions)
%         condition = Conditions{iC};
%         if isfield(normalized_mean_cycles, condition) && isfield(normalized_mean_cycles.(condition), muscle)
%             has_data = true;
%             break;
%         end
%     end
% 
%     if ~has_data
%         fprintf('Aucune donnée à afficher pour le muscle %s - graphique ignoré\n', muscle);
%         continue;
%     end
% 
%     % Créer une figure pour ce muscle
%     figure;
%     hold on;
% 
%     % Couleurs pour chaque condition
%     colors = {'b', 'g', 'r'};  % Bleu pour Plat, Vert pour Medium, Rouge pour High
%     legends = {};  % Pour stocker les légendes
%     h_legend = [];  % Pour stocker les handles des courbes de la légende
% 
%     % Parcourir chaque condition
%     for iC = 1:length(Conditions)
%         condition = Conditions{iC};
% 
%         if isfield(normalized_mean_cycles, condition) && isfield(normalized_mean_cycles.(condition), muscle)
%             % Récupérer le cycle normalisé pour cette condition
%             normalized_cycle = normalized_mean_cycles.(condition).(muscle);
% 
%             % Calculer l'écart type normalisé en gérant les NaN
%             if isfield(mean_cycles_combined, condition) && isfield(mean_cycles_combined.(condition), muscle)
%                 % Calculer la valeur max_flat pour la normalisation de l'écart type
%                 if isfield(global_mean_cycles, 'Plat') && isfield(global_mean_cycles.Plat, muscle)
%                     max_flat = max(global_mean_cycles.Plat.(muscle));
%                     if ~isnan(max_flat) && max_flat > 0
%                         std_cycle = nanstd(mean_cycles_combined.(condition).(muscle), 0, 1) / max_flat;
%                     else
%                         std_cycle = zeros(size(normalized_cycle));
%                     end
%                 else
%                     std_cycle = zeros(size(normalized_cycle));
%                 end
%             else
%                 std_cycle = zeros(size(normalized_cycle));
%             end
% 
%             % Tracer la courbe de la moyenne normalisée pour cette condition
%             h = plot(normalized_cycle, 'LineWidth', 2, 'Color', colors{iC}, 'DisplayName', condition);
% 
%             % Afficher la zone d'écart type autour de la moyenne (seulement si std_cycle est valide)
%             if ~any(isnan(std_cycle))
%                 fill([1:length(normalized_cycle), fliplr(1:length(normalized_cycle))], ...
%                     [normalized_cycle + std_cycle, fliplr(normalized_cycle - std_cycle)], ...
%                     colors{iC}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
%             end
% 
%             % Ajouter à la légende
%             legends{end+1} = [condition];
%             h_legend(end+1) = h;
%         end
%     end
% 
%     % Ajouter des légendes et un titre seulement s'il y a des données à afficher
%     if ~isempty(h_legend)
%         legend(h_legend, legends, 'Location', 'best');
%         title(['Cycle normalisé ' muscle ' - Toutes jambes - Population Adulte']);
%         xlabel('Cycle de marche (%)');
%         ylabel('Activation musculaire normalisée');
%         grid on;
% 
%         % % Ajouter une ligne horizontale à y=1 pour référence
%         % line([1, size(normalized_mean_cycles.Plat.(muscle), 2)], [1, 1], 'LineStyle', '--', 'Color', 'k');
% 
%         hold off;
% 
%         % Enregistrer la figure dans le dossier de résultats
%         saveas(gcf, fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ['Cycle_Normalisé_' muscle '_AllLegs.png']));
%         fprintf('Graphique sauvegardé pour le muscle %s\n', muscle);
%     else
%         close(gcf); % Fermer la figure vide
%         fprintf('Aucune courbe à afficher pour le muscle %s - figure fermée\n', muscle);
%     end
% end

% Créer une matrice structurée pour l'analyse de synergie
% Chaque colonne = un muscle, chaque ligne = un point du cycle
SYNERGY_Adulte = struct();

for iC = 1:length(Conditions)
    condition = Conditions{iC};
    
    % Initialiser la matrice pour cette condition
    condition_matrix = [];
    muscle_names = {};
    
    % Parcourir chaque muscle pour construire la matrice
    for iM = 1:length(Muscles)
        muscle = Muscles{iM};
        
        % Vérifier si ce muscle existe dans les données normalisées pour cette condition
        if isfield(normalized_mean_cycles, condition) && isfield(normalized_mean_cycles.(condition), muscle)
            % Récupérer le cycle normalisé (vecteur ligne)
            cycle_data = normalized_mean_cycles.(condition).(muscle);
            
            % Ajouter comme colonne à la matrice (transposer pour avoir en colonne)
            condition_matrix = [condition_matrix, cycle_data'];
            muscle_names{end+1} = muscle;
            
            fprintf('Muscle %s ajouté à la matrice pour la condition %s\n', muscle, condition);
        else
            fprintf('Muscle %s absent pour la condition %s - ignoré dans la matrice\n', muscle, condition);
        end
    end
    
    % Stocker la matrice et les noms des muscles pour cette condition
    if ~isempty(condition_matrix)
        SYNERGY_Adulte.(condition).matrix = condition_matrix;
        SYNERGY_Adulte.(condition).muscle_names = muscle_names;
        
        fprintf('Matrice créée pour %s : %d points x %d muscles\n', condition, size(condition_matrix, 1), size(condition_matrix, 2));
    else
        fprintf('Aucune donnée disponible pour créer la matrice de la condition %s\n', condition);
    end
end

suffix = upper(groupe_a_etudier);  % ex: 'ADULTES'

save(['SPM-EMG-' suffix '.mat'], 'mean_cycles_combined');
save(['SYNERGY-' suffix '.mat'], 'SYNERGY_Adulte');
save(['PARTICIPANT-CYCLES-' suffix '.mat'], 'normalized_participant_cycles', 'participant_tracking');

% Résumé de la structure créée
fprintf('\n=== RÉSUMÉ DE LA MATRICE POUR SYNERGIE ===\n');
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    if isfield(SYNERGY_Adulte, condition)
        fprintf('Condition %s:\n', condition);
        fprintf('  - Dimensions: %d points x %d muscles\n', ...
                size(SYNERGY_Adulte.(condition).matrix, 1), ...
                size(SYNERGY_Adulte.(condition).matrix, 2));
        fprintf('  - Muscles inclus: %s\n', strjoin(SYNERGY_Adulte.(condition).muscle_names, ', '));
    else
        fprintf('Condition %s: Aucune donnée\n', condition);
    end
end

fprintf('\nTraitement terminé. Vérifiez les messages ci-dessus pour identifier les données problématiques.\n');