function plot_all_bands(TABLES, pair_name, sous_phase, groupe_age)
    bands = {'Alpha', 'Beta', 'Gamma'};
    conditions = {'Plat', 'Medium', 'High'};
    
    % Couleurs pour chaque bande
    band_colors = [0.2 0.6 1;     % Alpha - Bleu
                   0.8 0.4 0.2;   % Beta - Orange foncé
                   0.6 0.2 0.8];  % Gamma - Violet
    
    figure('Position', [100 100 1200 400]);
    
    for iB = 1:length(bands)
        band = bands{iB};
        
        subplot(1, 3, iB);
        hold on;
        
        % Extraire les données pour cette bande
        T = TABLES.(pair_name).(band).(sous_phase);
        mask_groupe = ismember(T.GroupeAge, groupe_age);
        T_groupe = T(mask_groupe, :);
        
        % Tracer pour chaque condition
        for iC = 1:length(conditions)
            cond = conditions{iC};
            values = T_groupe.(cond);
            values = values(~isnan(values));
            
            if ~isempty(values)
                % Jitter pour séparer les points
                jitter = 0.15 * (rand(length(values), 1) - 0.5);
                x_pos = iC + jitter;
                
                % Points individuels
                scatter(x_pos, values, 50, band_colors(iB,:), 'filled', ...
                       'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', band_colors(iB,:)*0.7);
                
                % Moyenne
                mean_val = mean(values);
                line([iC-0.15, iC+0.15], [mean_val, mean_val], ...
                     'Color', band_colors(iB,:)*0.7, 'LineWidth', 2);
            end
        end
        
        % Mise en forme
        xlim([0.5 3.5]);
        xticks(1:3);
        xticklabels(conditions);
        title(sprintf('Bande %s', band), 'FontWeight', 'bold');
        
        if iB == 1
            ylabel('Cohérence EMG-EMG', 'FontWeight', 'bold');
        end
        if iB == 2
            xlabel('Conditions de surface', 'FontWeight', 'bold');
        end
        
        grid on;
        set(gca, 'GridAlpha', 0.3);
        
        % Ligne de référence
        yline(0.5, '--k', 'Alpha', 0.7);
    end
    
    sgtitle(sprintf('Cohérence %s - %s (%s)', ...
            strrep(pair_name,'_',' vs '), sous_phase, groupe_age), ...
            'FontSize', 14, 'FontWeight', 'bold');
end