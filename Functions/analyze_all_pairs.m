function analyze_all_pairs(TABLES, groupe_age, band)
    pairs = {'TAprox_TAdist', 'VL_RF', 'GM_SOL', 'GMED_RF', 'RF_ST'};
    
    for i = 1:length(pairs)
        fprintf('\n========================================\n');
        fprintf('ANALYSE DE LA PAIRE: %s\n', pairs{i});
        fprintf('========================================\n');
        
        % Modifier les paramètres dans le script principal
        assignin('base', 'pair_name', pairs{i});
        assignin('base', 'band', band);
        assignin('base', 'groupe_age', groupe_age);
        
        % Relancer l'analyse pour cette paire
        evalin('base', 'run(''script_principal'')');
        
        pause(2); % Pause pour voir les résultats
    end
end