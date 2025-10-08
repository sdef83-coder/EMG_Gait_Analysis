function generate_summary_report(TABLES, muscle_phases, groupe_age, bands)
    fprintf('\n=== RESUME STATISTIQUE ===\n');
    
    pair_names = fieldnames(muscle_phases);
    conditions = {'Plat', 'Medium', 'High'};
    
    % Créer un tableau de résumé
    summary_data = [];
    summary_labels = {};
    
    for iPair = 1:length(pair_names)
        pair_name = pair_names{iPair};
        sous_phases_interest = muscle_phases.(pair_name);
        
        for iP = 1:length(sous_phases_interest)
            sous_phase = sous_phases_interest{iP};
            
            for iBand = 1:length(bands)
                band = bands{iBand};
                
                % Vérifier si les données existent
                if ~isfield(TABLES.(pair_name), band) || ~isfield(TABLES.(pair_name).(band), sous_phase)
                    continue;
                end
                
                T = TABLES.(pair_name).(band).(sous_phase);
                mask_groupe = ismember(T.GroupeAge, groupe_age);
                T_groupe = T(mask_groupe, :);
                
                if height(T_groupe) == 0
                    continue;
                end
                
                % Calculer les moyennes pour chaque condition
                row_data = [];
                for iC = 1:length(conditions)
                    cond = conditions{iC};
                    values = T_groupe.(cond);
                    values = values(~isnan(values));
                    
                    if ~isempty(values)
                        row_data = [row_data, mean(values)];
                    else
                        row_data = [row_data, NaN];
                    end
                end
                
                summary_data = [summary_data; row_data];
                summary_labels{end+1} = sprintf('%s_%s_%s', pair_name, sous_phase, band);
            end
        end
    end
    
    % Afficher le tableau de résumé
    if ~isempty(summary_data)
        fprintf('\nTableau des moyennes de cohérence:\n');
        fprintf('%-25s | %-8s | %-8s | %-8s\n', 'Paire_Phase_Bande', 'Plat', 'Medium', 'High');
        fprintf('%s\n', repmat('-', 1, 55));
        
        for i = 1:size(summary_data, 1)
            fprintf('%-25s | %8.4f | %8.4f | %8.4f\n', ...
                   summary_labels{i}, summary_data(i,1), summary_data(i,2), summary_data(i,3));
        end
    end
end