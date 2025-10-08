clear;
close all;
clc;

addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\coherence_analysis'))
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Data\adults\');

[num, txt, raw] = xlsread('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\PhD projects\2) Projet_Surfaces_Irr\SCRIPTS\ActivationMusculaire\Mapping-EMG.xlsm');

Participant = {'CTL_25'}; % Exemple 1 seul participant
Condition = {'Plat',1;'Medium',3; 'High',5}; % Les numéros correspondent à leur décalage à droite dans le tableau de mapping .csv quant au 'CTL_??'
Essai = {'01' '02' '03' '04' '05' '06' '07' '08' '09' '10'};
raw_indices = struct( ...
    'TAprox', 13, ...
    'TAdist', 14, ...
    'SOL',    15, ...
    'GM',     16, ...
    'VL',     17, ...
    'RF',     18, ...
    'ST',     19 ...
);

alpha = 0.05;  % Seuil de significativité

% Paramètres de la transformation en ondelettes
FreqMin = 1;
FreqMax = 400;
Resolution = 1;
WaveNumber = 7;
FreqBands = struct( ...
    'Alpha', [8 12], ...
    'Beta', [13 30], ...
    'Gamma', [31 60] ...
);

% Nouvelle structure pour stocker les pourcentages de toe-off par condition/côté
ToeOffPercentages = struct();

% Extraction des données pour chaque pair de muscle, à chaque essai, pour chaque sujet
for iP = 1:length(Participant)
    load(['C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Matrix\' Participant{iP} '_MATRIX.mat'])
    run('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Association.m'); % Chemin vers script "Association.m"

    % Créer le dossier de sauvegarde
output_dir = fullfile('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Coherence', Participant{iP});
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

    Pairs = {}; % Réinitialise

    muscle_pairs = {
        'TAprox', 'TAdist';
        'VL', 'RF';
        'GM', 'SOL';
        'RF', 'ST';
        % Ajoute ici d'autres paires si nécessaire
    };

    for i = 1:size(muscle_pairs, 1)
        m1 = muscle_pairs{i,1};
        m2 = muscle_pairs{i,2};

        for side = {'left', 'right'}
            side_str = side{1};

            if strcmp(side_str, 'left')
                assoc = sensor_association_left;
            else
                assoc = sensor_association_right;
            end

            % Vérifie que les deux muscles sont bien présents dans la structure
            if isfield(assoc, ['EMG_' m1]) && isfield(assoc, ['EMG_' m2]) ...
                    && isfield(raw_indices, m1) && isfield(raw_indices, m2)

                % Ajoute à la cellule Pairs : sensor1, sensor2, m1, m2, côté, ligne m1, ligne m2
                Pairs(end+1,:) = {
                    assoc.(['EMG_' m1]), ...
                    assoc.(['EMG_' m2]), ...
                    m1, ...
                    m2, ...
                    side_str, ...
                    raw_indices.(m1), ...
                    raw_indices.(m2)
                };
            end
        end
    end

    % CALCUL DES POURCENTAGES MOYENS DE TOE-OFF PAR CONDITION/CÔTÉ
    for iC = 1:size(Condition,1)
        for side_temp = {'left', 'right'}
            side_str_temp = side_temp{1};
            
            if strcmp(side_str_temp, 'left')
                Side_temp = {'left',0} ;
            else
                Side_temp = {'right',1} ;
            end
            
            % Initialiser les variables pour calculer le % moyen de toe-off
            toeoff_percentages = [];
            
            for iEs = 1:length(Essai)
                file = [Participant{iP} '_' Condition{iC,1} '_' Essai{iEs} '.c3d'];
                if ~isfile(file)
                    continue;
                end
                
                data = btkReadAcquisition(file);
                analogs = btkGetAnalogs(data);
                FreqS = btkGetAnalogFrequency(data);
                FreqVicon = 100;
                
                % Obtenir un signal EMG temporaire pour les fonctions de détection
                fieldnames_analogs = fieldnames(analogs);
                temp_emg = analogs.(fieldnames_analogs{1});
                EMG_temp = [temp_emg, temp_emg]; % Format attendu par vos fonctions
                
                % Obtenir les indices selon le côté
                if strcmp(side_str_temp, 'left')
                    heel_strikes = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_temp);
                    toe_offs = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                else
                    heel_strikes = indiceRight(data, analogs, FreqS, FreqVicon, EMG_temp);
                    toe_offs = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_temp); % Supposée fonction équivalente
                end
                
                % Calculer les pourcentages de toe-off pour chaque cycle
                for i = 1:min(length(heel_strikes)-1, length(toe_offs))
                    if toe_offs(i) > heel_strikes(i) && toe_offs(i) < heel_strikes(i+1)
                        cycle_length = heel_strikes(i+1) - heel_strikes(i);
                        toeoff_relative = toe_offs(i) - heel_strikes(i);
                        toeoff_percent = (toeoff_relative / cycle_length) * 100;
                        if toeoff_percent > 0 && toeoff_percent < 100 % Validation
                            toeoff_percentages(end+1) = toeoff_percent;
                        end
                    end
                end
            end
            
            % Calculer le pourcentage moyen de toe-off pour cette condition/côté
            if ~isempty(toeoff_percentages)
                ToeOffPercentages.(Condition{iC,1}).(Side_temp{1,1}) = mean(toeoff_percentages);
                fprintf('Toe-off moyen calculé pour %s - %s: %.1f%%\n', Condition{iC,1}, Side_temp{1,1}, mean(toeoff_percentages));
            else
                ToeOffPercentages.(Condition{iC,1}).(Side_temp{1,1}) = 60; % Valeur par défaut
                fprintf('Toe-off par défaut utilisé pour %s - %s: 60%%\n', Condition{iC,1}, Side_temp{1,1});
            end
        end
    end

    % TRAITEMENT PRINCIPAL
    for iC = 1:size(Condition,1)
        for iPair = 1:size(Pairs, 1)
            if strcmp(Pairs{iPair,5},'left')
                Side = {'left',0} ;
            else
                Side = {'right',1} ;
            end

            %% Verifie si un des deux muscles est a rejeter
            for ii = 1:length(raw)
                if strcmp(raw{1,ii},Participant{iP})
                    break
                end
            end

            if raw{Pairs{iPair,6},ii+Condition{iC,2}+Side{1,2}}==2 &&... 
                    raw{Pairs{iPair,7},ii+Condition{iC,2}+Side{1,2}}==2

                %% Identifie quels sont les cycles a supprimer 
                BadCycles = CYCLES_OUTLIERS.(Participant{iP}).(Condition{iC,1}).(Side{1,1});
                tutu = BadCycles.(['EMG_' Pairs{iPair,3}]) + BadCycles.(['EMG_' Pairs{iPair,4}]);
         
                %% Initialiser les matrices qui vont stocker les moyennes
                PowSpec_s1 = zeros(400,1000);
                PowSpec_s2 = zeros(400,1000);
                cross_spectrum = zeros(400,1000);
                NcycleP = 0; NcycleC = 0;

                for iEs = 1:length(Essai)
                    % Charger le fichier pour l'essai actuel
                    file = [Participant{iP} '_' Condition{iC,1} '_' Essai{iEs} '.c3d'];
                    if ~isfile(file)
                        warning(['Fichier non trouvé : ' file]);
                        continue;
                    end

                    % Lire les données c3d
                    data = btkReadAcquisition(file);
                    analogs = btkGetAnalogs(data);

                    % Réinitialiser les matrices à chaque paire de muscles
                    EMG = [];
                    EMG(:, 1) = analogs.(Pairs{iPair, 1});
                    EMG(:, 2) = analogs.(Pairs{iPair, 2});
                    muscle1_name = Pairs{iPair, 3};
                    muscle2_name = Pairs{iPair, 4};
                    FreqS = btkGetAnalogFrequency(data);
                    FreqVicon = 100;

                    % Traiter les cycles
                    if strcmp(Side{1}, 'left')
                        Cycles = indiceLeft(data, analogs, FreqS, FreqVicon, EMG);
                    else
                        Cycles = indiceRight(data, analogs, FreqS, FreqVicon, EMG);
                    end
                    Args = WaveletParameters(FreqMin, FreqMax, Resolution, WaveNumber, FreqS);

                    % Filtrer le signal EMG
                    EMG = filtrage(EMG, FreqS, 8, 400);
                    EMG = EMG - repmat(mean(EMG), length(EMG), 1);

                    % Calcul de la carte temps-fréquence - Premier calcul pour obtenir les dimensions correctes
                    [TFR(:,:,1), period, ~, ~] = wavelet(EMG(:,1), Args.DT, Args.Pad, Args.DJ, Args.S0, Args.J1, Args.Mother, Args.Cycles);
                    [TFR(:,:,2), ~, ~, ~]      = wavelet(EMG(:,2), Args.DT, Args.Pad, Args.DJ, Args.S0, Args.J1, Args.Mother, Args.Cycles);
                    Freq = 1 ./ period; % Mettre à jour Freq avec les périodes réelles

                    % Découpage cycle par cycle 
                    for iCycles = 1:length(Cycles) - 1
                        NcycleP = NcycleP + 1;
                        if tutu(NcycleP)==0  % VOTRE CONDITION DE FILTRAGE 

                            TFR_cycle(:,:,1) = TFR(:, Cycles(iCycles):Cycles(iCycles+1), 1); % Découper les cycles
                            TFR_cycle(:,:,2) = TFR(:, Cycles(iCycles):Cycles(iCycles+1), 2); % Découper les cycles

                            [X, Y] = meshgrid(1:size(TFR_cycle, 2), 1:size(TFR_cycle, 1)); % Interpolation des cycles temps-fréquence
                            [Xq, Yq] = meshgrid(linspace(1, size(TFR_cycle, 2), 1000), 1:size(TFR_cycle, 1));
                            TFR_int(:,:,1) = interp2(X, Y, TFR_cycle(:,:,1), Xq, Yq, 'spline');
                            TFR_int(:,:,2) = interp2(X, Y, TFR_cycle(:,:,2), Xq, Yq, 'spline');

                            PowSpec_s1 = PowSpec_s1 + abs(TFR_int(:,:,1)).^2;
                            PowSpec_s2 = PowSpec_s2 + abs(TFR_int(:,:,2)).^2;
                            cross_spectrum = cross_spectrum + (TFR_int(:,:,1)) .* conj(TFR_int(:,:,2));

                            NcycleC = NcycleC + 1;
                            clear TFR_int TFR_cycle
                        end
                    end
                    clear TFR
                end

                % Moyenne et cohérence
                if NcycleC == 0
                    warning(['Aucun cycle retenu pour ' muscle1_name ' - ' muscle2_name ' | ' ...
                        Condition{iC,1} ' | ' Side{1,1}]);
                    continue;
                end

                PowSpec_s1 = PowSpec_s1 / NcycleC;
                PowSpec_s2 = PowSpec_s2 / NcycleC;
                cross_spectrum_mean = cross_spectrum / NcycleC;
                cross_spectrum_abs = abs(cross_spectrum_mean).^2;
                Coherence = cross_spectrum_abs ./ (PowSpec_s1 .* PowSpec_s2);

                % SEGMENTATION PAR PHASES
                % Obtenir le pourcentage de toe-off pour cette condition/côté
                toeoff_percent = ToeOffPercentages.(Condition{iC,1}).(Side{1,1});
                toeoff_index = round(toeoff_percent * 10); % Sur 1000 points, toeoff_percent * 10

                % S'assurer que l'index est dans les limites
                toeoff_index = max(1, min(toeoff_index, 999));

                % Définir les phases
                stance_phase_indices = 1:toeoff_index;                    % Phase d'appui
                swing_phase_indices = (toeoff_index+1):1000;              % Phase d'oscillation

                % Calcul de la cohérence moyenne par bande de fréquence et par phase
                mean_coherence_stance = struct();
                mean_coherence_swing = struct();
                mean_coherence_full = struct();

                for bandName = fieldnames(FreqBands)'
                    band = bandName{1};
                    range = FreqBands.(band);
                    idx_band = find(Freq >= range(1) & Freq <= range(2));
                    
                    % Cohérence pour le cycle complet (COMME AVANT)
                    mean_coherence_full.(band) = mean(mean(Coherence(idx_band, :), 2), 'omitnan');
                    
                    % Cohérence pour la phase d'appui
                    mean_coherence_stance.(band) = mean(mean(Coherence(idx_band, stance_phase_indices), 2), 'omitnan');
                    
                    % Cohérence pour la phase d'oscillation
                    mean_coherence_swing.(band) = mean(mean(Coherence(idx_band, swing_phase_indices), 2), 'omitnan');
                    
                    % Sauvegarder dans DATA
                    DATA.(Condition{iC,1}).(Side{1,1}).(['MeanCoherence_' band '_' muscle1_name '_' muscle2_name]) = mean_coherence_full.(band); % VOTRE FORMAT CONSERVÉ
                    DATA.(Condition{iC,1}).(Side{1,1}).(['MeanCoherence_Stance_' band '_' muscle1_name '_' muscle2_name]) = mean_coherence_stance.(band);
                    DATA.(Condition{iC,1}).(Side{1,1}).(['MeanCoherence_Swing_' band '_' muscle1_name '_' muscle2_name]) = mean_coherence_swing.(band);
                end

                % Sauvegarder le pourcentage de toe-off utilisé
                DATA.(Condition{iC,1}).(Side{1,1}).(['ToeOffPercent_' muscle1_name '_' muscle2_name]) = toeoff_percent;

                seuil = 1 - 0.05^(1 / (NcycleC - 1));

                % Visualisation améliorée avec les phases
                Time = linspace(0, 100, 1000);
                figure;
                
                subplot(2, 3, 1); 
                imagesc(Time, Freq, PowSpec_s1); 
                title(['Power Spec ' muscle1_name]);
                xlabel('% Cycle'); ylabel('Fréquence (Hz)');
                
                subplot(2, 3, 2); 
                imagesc(Time, Freq, PowSpec_s2); 
                title(['Power Spec ' muscle2_name]);
                xlabel('% Cycle'); ylabel('Fréquence (Hz)');
                
                subplot(2, 3, 3); 
                imagesc(Time, Freq, Coherence > seuil); 
                title([muscle1_name '-' muscle2_name ' Coherence']);
                hold on;
                % Ligne verticale pour marquer le toe-off
                line([toeoff_percent toeoff_percent], [min(Freq) max(Freq)], 'Color', 'red', 'LineWidth', 2);
                xlabel('% Cycle'); ylabel('Fréquence (Hz)');
                
                subplot(2, 3, 4); 
                plot(Freq, mean(Coherence, 2)); 
                hold on;
                line([0 400], [seuil seuil], 'Color', 'red');
                xlim([1 400]); ylim([0 1]);
                title('Cohérence - Cycle complet');
                xlabel('Fréquence (Hz)'); ylabel('Cohérence');
                
                subplot(2, 3, 5); 
                plot(Freq, mean(Coherence(:, stance_phase_indices), 2)); 
                hold on;
                line([0 400], [seuil seuil], 'Color', 'red');
                xlim([1 400]); ylim([0 1]);
                title('Cohérence - Phase d''appui');
                xlabel('Fréquence (Hz)'); ylabel('Cohérence');
                
                subplot(2, 3, 6); 
                plot(Freq, mean(Coherence(:, swing_phase_indices), 2)); 
                hold on;
                line([0 400], [seuil seuil], 'Color', 'red');
                xlim([1 400]); ylim([0 1]);
                title('Cohérence - Phase oscillation');
                xlabel('Fréquence (Hz)'); ylabel('Cohérence');
                
                sgtitle(['Condition : ' Condition{iC,1} ' | Côté : ' Side{1,1} ' | ' muscle1_name ' - ' muscle2_name ' | Toe-off: ' num2str(toeoff_percent, '%.1f') '%']);

                fig_name = sprintf('%s_%s_%s-%s', Condition{iC,1}, Side{1,1}, muscle1_name, muscle2_name);
                saveas(gcf, fullfile(output_dir, [fig_name '.png']));
                close(gcf);

                % Affichage des résultats dans la console
                fprintf('\n=== Résultats pour %s-%s | %s | %s ===\n', muscle1_name, muscle2_name, Condition{iC,1}, Side{1,1});
                fprintf('Toe-off moyen : %.1f%% du cycle\n', toeoff_percent);
                fprintf('Phase d''appui : 0-%.1f%% | Phase oscillation : %.1f-100%%\n', toeoff_percent, toeoff_percent);
                
                for bandName = fieldnames(FreqBands)'
                    band = bandName{1};
                    fprintf('%s - Full: %.3f | Stance: %.3f | Swing: %.3f\n', ...
                        band, mean_coherence_full.(band), mean_coherence_stance.(band), mean_coherence_swing.(band));
                end

                % Sauvegarde des résultats dans la structure
                DATA.(Condition{iC,1}).(Side{1,1}).(['PowSpec' muscle1_name]) = PowSpec_s1;
                DATA.(Condition{iC,1}).(Side{1,1}).(['PowSpec' muscle2_name]) = PowSpec_s2;
                DATA.(Condition{iC,1}).(Side{1,1}).(['Cross_Spectrum_' muscle1_name '_' muscle2_name]) = cross_spectrum_abs;
                DATA.(Condition{iC,1}).(Side{1,1}).(['Coherence_' muscle1_name '_' muscle2_name]) = Coherence;
                DATA.(Condition{iC,1}).(Side{1,1}).(['Ncycle_' muscle1_name '_' muscle2_name]) = NcycleC;
                DATA.(Condition{iC,1}).(Side{1,1}).(['Seuil_' muscle1_name '_' muscle2_name])  = seuil;
                DATA.(Condition{iC,1}).(Side{1,1}).Freq = Freq;
                
                % NOUVELLES SAUVEGARDES POUR LES PHASES
                DATA.(Condition{iC,1}).(Side{1,1}).(['Coherence_Stance_' muscle1_name '_' muscle2_name]) = Coherence(:, stance_phase_indices);
                DATA.(Condition{iC,1}).(Side{1,1}).(['Coherence_Swing_' muscle1_name '_' muscle2_name]) = Coherence(:, swing_phase_indices);
                DATA.(Condition{iC,1}).(Side{1,1}).(['StanceIndices_' muscle1_name '_' muscle2_name]) = stance_phase_indices;
                DATA.(Condition{iC,1}).(Side{1,1}).(['SwingIndices_' muscle1_name '_' muscle2_name]) = swing_phase_indices;
            end
        end
    end
end

save(fullfile(output_dir, ['Coherence_' Participant{iP} '.mat']), 'DATA');