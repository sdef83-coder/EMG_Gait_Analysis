function plot_single_phase(T_groupe, conditions, condition_colors, sous_phase, pair_name, band, groupe_age, is_subplot)
    if nargin < 8
        is_subplot = false;
    end

    hold on;
    
    % Paramètres de visualisation
    jitter_width = 0.15;
    point_size = 60;
    alpha_value = 0.7;
    
    stats_data = [];
    
    for iC = 1:length(conditions)
        cond = conditions{iC};
        
        % Extraire les valeurs pour cette condition
        values = T_groupe.(cond);
        values = values(~isnan(values));
        
        if ~isempty(values)
            % Position X de base pour cette condition
            x_base = iC;
            
            % Ajouter du jitter
            jitter = jitter_width * (rand(length(values), 1) - 0.5);
            x_pos = x_base + jitter;
            
            % Tracer les points individuels
            scatter(x_pos, values, point_size, condition_colors(iC,:), 'filled', ...
                   'MarkerFaceAlpha', alpha_value, 'MarkerEdgeColor', condition_colors(iC,:)*0.7, ...
                   'LineWidth', 1);
            
            % Calculer et afficher la moyenne
            mean_val = mean(values);
            std_val = std(values);
            
            % Ligne horizontale pour la moyenne
            line([x_base-0.2, x_base+0.2], [mean_val, mean_val], ...
                 'Color', condition_colors(iC,:)*0.7, 'LineWidth', 3);
            
            % Barres d'erreur (écart-type)
            errorbar(x_base, mean_val, std_val, 'Color', condition_colors(iC,:)*0.5, ...
                    'LineWidth', 2, 'CapSize', 10, 'LineStyle', 'none');
            
            fprintf('%s: n=%d, moyenne=%.3f ± %.3f\n', cond, length(values), mean_val, std_val);
            stats_data = [stats_data; values];
        else
            fprintf('%s: Aucune donnée disponible\n', cond);
        end
    end
    
    % Mise en forme
    xlim([0.5 3.5]);
    
    if ~isempty(stats_data)
        y_min = min(stats_data) - 0.1;
        y_max = max(stats_data) + 0.2;
        ylim([y_min y_max]);
        
        % Ligne de référence (Décommenter si besoin de rajouter un seuil)
        % yline(0.5, '--k', 'LineWidth', 1.5, 'Alpha', 0.7,'Label', 'Seuil significatif', 'LabelHorizontalAlignment', 'right');
    end
    
    xticks(1:3);
    xticklabels(conditions);
    
    if ~is_subplot
        xlabel('Conditions de surface', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel('Cohérence EMG-EMG', 'FontSize', 12, 'FontWeight', 'bold');
        title(sprintf('Cohérence %s : %s - %s (%s, n=%d)', ...
              band, strrep(pair_name,'_',' vs '), sous_phase, groupe_age, height(T_groupe)), ...
              'FontSize', 14, 'FontWeight', 'bold');
    end
    
    grid on;
    set(gca, 'GridAlpha', 0.3, 'FontSize', 11, 'LineWidth', 1);
    box on;
end