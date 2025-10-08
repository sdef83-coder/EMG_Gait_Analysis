function plot_all_bands_specific(TABLES, pair_name, groupe_age)
    % Configuration spécifique par paire de muscles
    muscle_phases = struct();
    muscle_phases.TAprox_TAdist = {'LoadingResponse', 'PreSwing'};
    muscle_phases.VL_RF = {'LoadingResponse'};
    muscle_phases.GM_SOL = {'MidStance'};
    muscle_phases.GMED_RF = {'LoadingResponse'};
    muscle_phases.RF_ST = {'LoadingResponse'};
    
    if ~isfield(muscle_phases, pair_name)
        error('Paire de muscles "%s" non définie', pair_name);
    end
    
    sous_phases_interest = muscle_phases.(pair_name);
    bands = {'Alpha', 'Beta', 'Gamma'};
    conditions = {'Plat', 'Medium', 'High'};
    
    % Couleurs pour chaque bande
    band_colors = [0.2 0.6 1;     % Alpha - Bleu
                   0.8 0.4 0.2;   % Beta - Orange foncé
                   0.6 0.2 0.8];  % Gamma - Violet
    
    num_phases = length(sous_phases_interest);
    
    % Créer une figure adaptée au nombre de sous-phases et bandes
    figure('Position', [100 100 400*length(bands) 300*num_phases]);
    
    subplot_idx = 1;
    
    for iP = 1:num_phases
        sous_phase = sous_phases_interest{iP};
        
        for iB = 1:length(bands)
            band = bands{iB};
            
            subplot(num_phases, length(bands), subplot_idx);
            hold on;
            
            % Extraire les données pour cette bande et sous-phase
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
            
            % Titre avec bande et sous-phase
            if num_phases == 1
                title(sprintf('%s', band), 'FontWeight', 'bold');
            else
                title(sprintf('%s - %s', band, sous_phase), 'FontWeight', 'bold');
            end
            
            % Labels des axes
            if iB == 1 && iP == 1
                ylabel('Cohérence EMG-EMG', 'FontWeight', 'bold');
            end
            if iP == num_phases && iB == 2
                xlabel('Conditions de surface', 'FontWeight', 'bold');
            end
            
            grid on;
            set(gca, 'GridAlpha', 0.3);
            yline(0.5, '--k', 'Alpha', 0.7);
            
            subplot_idx = subplot_idx + 1;
        end
    end
    
    sgtitle(sprintf('Cohérence %s (%s) - Toutes bandes', ...
            strrep(pair_name,'_',' vs '), groupe_age), ...
            'FontSize', 14, 'FontWeight', 'bold');
end