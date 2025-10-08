% Fonction complémentaire pour caractériser les artefacts
function [artifacts_info] = characterizeArtifacts(signal_filtered, sampling_rate, varargin)
    % CHARACTERIZEARTIFACTS analyse les caractéristiques des artefacts dans un signal_filtered
    %
    % Syntaxe :
    %   artifacts_info = characterizeArtifacts(signal_filtered, sampling_rate)
    %   artifacts_info = characterizeArtifacts(signal_filtered, sampling_rate, 'Param1', Valeur1, ...)
    %
    % Entrées :
    %   signal_filtered        - Vecteur du signal_filtered à analyser
    %   sampling_rate - Fréquence d'échantillonnage en Hz
    %
    % Paramètres optionnels (paires 'nom', valeur) :
    %   'Method'           - Méthode de détection : 'MAD' (défaut), 'STD', 'IQR', 'Custom'
    %   'ThresholdFactor'  - Facteur multiplicatif pour le seuil (défaut: 3.5)
    %   'CustomThresholds' - [seuil_bas, seuil_haut] pour méthode 'Custom'
    %   'MinDuration'      - Durée minimale (en ms) pour considérer un artefact (défaut: 0)
    %   'MaxDuration'      - Durée maximale (en ms) pour considérer un artefact (défaut: Inf)
    %   'Visualize'        - Afficher les visualisations (true/false, défaut: true)
    %   'BinWidth'         - Largeur des bins pour l'histogramme en ms (défaut: auto)

    % Traiter les arguments d'entrée
    p = inputParser;
    addRequired(p, 'signal_filtered', @isnumeric);
    addRequired(p, 'sampling_rate', @isnumeric);
    addParameter(p, 'Method', 'STD', @ischar);
    addParameter(p, 'ThresholdFactor', 5, @isnumeric);
    addParameter(p, 'CustomThresholds', [], @isnumeric);
    addParameter(p, 'MinDuration', 0, @isnumeric);
    addParameter(p, 'MaxDuration', Inf, @isnumeric);
    addParameter(p, 'Visualize', true, @islogical);
    addParameter(p, 'BinWidth', [], @isnumeric);
    parse(p, signal_filtered, sampling_rate, varargin{:});

    % Extraire les paramètres
    method = p.Results.Method;
    threshold_factor = p.Results.ThresholdFactor;
    custom_thresholds = p.Results.CustomThresholds;
    min_duration_ms = p.Results.MinDuration;
    max_duration_ms = p.Results.MaxDuration;
    visualize = p.Results.Visualize;
    bin_width = p.Results.BinWidth;

    % Détection des artefacts selon la méthode choisie
    switch upper(method)
        case 'MAD'
            med_val = median(signal_filtered);
            mad_val = median(abs(signal_filtered - med_val));
            upper_threshold = med_val + threshold_factor * mad_val;
            lower_threshold = med_val - threshold_factor * mad_val;
            
        case 'STD'
            mean_val = mean(signal_filtered);
            std_val = std(signal_filtered);
            upper_threshold = mean_val + threshold_factor * std_val;
            lower_threshold = mean_val - threshold_factor * std_val;
            
        case 'IQR'
            q1 = prctile(signal_filtered, 25);
            q3 = prctile(signal_filtered, 75);
            iqr_val = q3 - q1;
            upper_threshold = q3 + threshold_factor * iqr_val;
            lower_threshold = q1 - threshold_factor * iqr_val;
            
        case 'CUSTOM'
            if length(custom_thresholds) ~= 2
                error('Pour la méthode CUSTOM, CustomThresholds doit être [seuil_bas, seuil_haut]');
            end
            lower_threshold = custom_thresholds(1);
            upper_threshold = custom_thresholds(2);
            
        otherwise
            error('Méthode de détection non reconnue');
    end

    % Détection des artefacts
    artifacts = (signal_filtered > upper_threshold) | (signal_filtered < lower_threshold);

    % Trouver les indices de début et fin des artefacts
    artifact_starts = find(diff([0; artifacts]) == 1);
    artifact_ends = find(diff([artifacts; 0]) == -1);

    % Calculer les durées
    if ~isempty(artifact_starts)
        artifact_durations_samples = artifact_ends - artifact_starts + 1;
        artifact_durations_ms = artifact_durations_samples / sampling_rate * 1000;
        
        % Filtrer par durée
        valid_artifacts = (artifact_durations_ms >= min_duration_ms) & ...
                          (artifact_durations_ms <= max_duration_ms);
        
        artifact_starts = artifact_starts(valid_artifacts);
        artifact_ends = artifact_ends(valid_artifacts);
        artifact_durations_samples = artifact_durations_samples(valid_artifacts);
        artifact_durations_ms = artifact_durations_ms(valid_artifacts);
    end

    % Statistiques des artefacts
    num_artifacts = length(artifact_starts);
    artifacts_info = struct();
    artifacts_info.method = method;
    artifacts_info.threshold_factor = threshold_factor;
    artifacts_info.lower_threshold = lower_threshold;
    artifacts_info.upper_threshold = upper_threshold;
    artifacts_info.count = num_artifacts;
    artifacts_info.indices = [artifact_starts, artifact_ends];
    artifacts_info.percentage = sum(artifacts)/length(artifacts)*100;

    if num_artifacts > 0
        artifacts_info.mean_duration_ms = mean(artifact_durations_ms);
        artifacts_info.median_duration_ms = median(artifact_durations_ms);
        artifacts_info.std_duration_ms = std(artifact_durations_ms);
        artifacts_info.max_duration_ms = max(artifact_durations_ms);
        artifacts_info.min_duration_ms = min(artifact_durations_ms);
        artifacts_info.durations_ms = artifact_durations_ms;
    else
        artifacts_info.mean_duration_ms = NaN;
        artifacts_info.median_duration_ms = NaN;
        artifacts_info.std_duration_ms = NaN;
        artifacts_info.max_duration_ms = NaN;
        artifacts_info.min_duration_ms = NaN;
        artifacts_info.durations_ms = [];
    end
% 
%     % Visualisations
%     if visualize && ~isempty(signal_filtered)
%         % Figure 1: signal_filtered avec artefacts
%         figure('Name', 'Détection des artefacts');
% 
%         % Subplot supérieur: signal_filtered et seuils
%         subplot(2,1,1);
%         t = (0:length(signal_filtered)-1) / sampling_rate;
%         plot(t, signal_filtered, 'b');
%         hold on;
%         plot([t(1), t(end)], [upper_threshold, upper_threshold], 'r--');
%         plot([t(1), t(end)], [lower_threshold, lower_threshold], 'r--');
% 
%         % Surligner les artefacts
%         for i = 1:num_artifacts
%             x_start = (artifact_starts(i)-1) / sampling_rate;
%             x_end = (artifact_ends(i)-1) / sampling_rate;
%             h = patch([x_start x_end x_end x_start], ...
%                   [min(signal_filtered) min(signal_filtered) max(signal_filtered) max(signal_filtered)], ...
%                   'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
%         end
% 
%         title('signal_filtered avec artefacts détectés');
%         xlabel('Temps (s)');
%         ylabel('Amplitude');
%         legend('signal_filtered', 'Seuils');
% 
%         % Subplot inférieur: histogramme des durées
%         subplot(2,1,2);
%         if num_artifacts > 0
%             if isempty(bin_width)
%                 histogram(artifact_durations_ms);
%             else
%                 histogram(artifact_durations_ms, 'BinWidth', bin_width);
%             end
%             title('Distribution de la durée des artefacts');
%             xlabel('Durée (ms)');
%             ylabel('Fréquence');
%         else
%             text(0.5, 0.5, 'Aucun artefact détecté', 'HorizontalAlignment', 'center');
%             set(gca, 'XTick', [], 'YTick', []);
%         end
%     end
% end