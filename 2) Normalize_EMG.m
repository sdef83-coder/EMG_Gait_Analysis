%% === EXTRAIRE LA MATRICE DE SYNERGY ET LES MATRICES INDIVIDUELLES ===
% Les matrices indiv. contiennent : 
% 1) Matrices de synergies (SYNERGY_MATRIX) ; 
% 2) Signaux de chaque cycle normalisés (CYCLES_SIGNAL_NORMALIZED) ; 
% 3) Cycles moyens bruts par muscle/jambe/condition (CYCLES_MOYENS_BRUTS) ;
% 4) Les cycles outliers (CYCLES_OUTLIERS) ;
% 5) Le nombre de cycles valides par condition/jambe (CYCLES_COUNT) ;
% 6) Les % du cycle des Toe-offs (CYCLES_TOEOFF)
%
% Normalisation CEDE-compatible sans MVC :
%  -> Référence unique "dans la tâche" (Plat) = ref_max_plat (max de tous les points
%     de tous les cycles valides sur Plat) calculée par (participant×jambe×muscle).
%  -> Même ref utilisée pour normaliser Plat/Medium/High, et pour SYNERGY_MATRIX et CYCLES_SIGNAL_NORMALIZED.

clc; clear; close all;

addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire'));

cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Data\jeunes_enfants\');

Participant = {'CTL_63'};
Condition   = {'Plat','Medium','High'};
Essai       = {'01','02','03','04'}; % '05' '06' ...
muscles     = {'EMG_TAprox','EMG_TAdist','EMG_SOL','EMG_GM','EMG_VL','EMG_RF','EMG_ST','EMG_GMED'};
jambes      = {'left','right'};

run Association.m

% === Structures de sortie ===
CYCLES_MOYENS_BRUTS     = struct();
CYCLES_MOYENS           = struct();
CYCLES_STD              = struct();
CYCLES_SIGNAL           = struct();
CYCLES_SIGNAL_NORMALIZED= struct();
CYCLES_OUTLIERS         = struct();
CYCLES_COUNT            = struct();
CYCLES_TOEOFF           = struct();
SYNERGY_MATRIX          = struct();

% Références de normalisation (unique par participant×jambe×muscle)
REF_MAX_PLAT            = struct();

% Désactiver l'affichage des figures pendant le traitement
original_visible = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

fprintf('=== DÉBUT DU TRAITEMENT ===\n');

for j = 1:length(jambes)
    fprintf('Traitement jambe: %s\n', jambes{j});

    % Association capteurs
    if strcmp(jambes{j}, 'left')
        sensor_association = sensor_association_left;
    else
        sensor_association = sensor_association_right;
    end

    for iP = 1:length(Participant)
        pid = Participant{iP};
        fprintf('  Participant: %s\n', pid);

        % === ÉTAPE 1 — TRAITER "PLAT" pour construire les références ===
        fprintf('  === ÉTAPE 1: Traitement de PLAT (référence unique) ===\n');
        iC_plat = find(strcmp(Condition,'Plat'),1);
        if isempty(iC_plat)
            error('Condition "Plat" introuvable dans la liste Condition.');
        end
        cond_ref = Condition{iC_plat};

        % Dossier
        participant_folder_plat = fullfile( ...
            'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ...
            pid, cond_ref);
        if ~exist(participant_folder_plat,'dir'); mkdir(participant_folder_plat); end

        fprintf('    Condition: %s (référence)\n', cond_ref);
        fprintf('      Extraction des cycles Plat + Toe-Off...\n');

        % On remplit CYCLES_SIGNAL.(pid).(Plat).(jambe).(muscle) = {cycles bruts interpolés}
        % + on calcule outliers par muscle
        for m = 1:length(muscles)
            muscle = muscles{m};

            figure_concatenated = figure('Name', sprintf('Cycles Concaténés - %s - %s - %s - %s', pid, cond_ref, jambes{j}, muscle));
            hold on;

            cycle_offset   = 0;
            cycle_positions= [];
            muscle_maxima  = [];
            total_cycles   = 0;
            all_cycles_data= {}; % Cellule des cycles (100 pts) bruts (enveloppe nettoyée)

            all_toeoff_percentages = []; % calcul une fois par essai au m==1 ci-dessous

            for iEs = 1:length(Essai)
                file = [pid '_' cond_ref '_' Essai{iEs} '.c3d'];
                try
                    data   = btkReadAcquisition(file);
                    analogs= btkGetAnalogs(data);

                    sensor_name = sensor_association.(muscle);
                    EMG_signal  = analogs.(sensor_name);

                    % Traitement signal
                    FreqS      = btkGetAnalogFrequency(data);
                    FreqVicon  = 100;
                    EMGfilt    = filtrage(EMG_signal, FreqS, 20, 450);
                    artifacts_info = characterizeArtifacts(EMGfilt, FreqS);
                    [signal_cleaned, EMGenvcleaned] = cleanEMGdata(EMGfilt, FreqS, 'STD', 5, 15, false); % Defini le low-pass à 10 dans cette fonction

                    % Événements
                    if strcmp(jambes{j}, 'left')
                        HS = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_signal);
                        TO = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                    else
                        HS = indiceRight(data, analogs, FreqS, FreqVicon, EMG_signal);
                        TO = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_signal);
                    end

                    num_cycles = length(HS)-1;
                    cycles_frames = cell(num_cycles,1);

                    % Toe-off % (calculer une seule fois par essai lorsque m==1)
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
                                warning('Pas de toe-off (Plat) : cycle %d, essai %s', ii, Essai{iEs});
                                all_toeoff_percentages = [all_toeoff_percentages, 60];
                            end
                        end
                    end

                    for ii = 1:num_cycles
                        cycles_frames{ii} = HS(ii):HS(ii+1);
                    end

                    for ii = 1:num_cycles
                        cdata = EMGenvcleaned(cycles_frames{ii});
                        cinterp = interp1(linspace(1,length(cdata),length(cdata)), cdata, linspace(1,length(cdata),100), 'pchip');

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

            % traits séparateurs
            y_limits = ylim;
            for ii = 1:length(cycle_positions)
                plot([cycle_positions(ii), cycle_positions(ii)], y_limits, 'r', 'LineWidth', 0.5);
            end
            title(sprintf('All cycles for %s - %s - Muscle: %s - Leg: %s', pid, cond_ref, muscle, jambes{j}), 'Interpreter','none');
            xlabel('Temps normalisé'); ylabel('Enveloppe');

            % Outliers (par muscle)
            outlier_indices = detect_outliers(muscle_maxima, cycle_offset, 3);

            % Save figure
            save_path = fullfile(participant_folder_plat, sprintf('%s_%s_%s_%s.png', pid, cond_ref, jambes{j}, muscle));
            print(figure_concatenated, save_path, '-dpng', '-r150'); close(figure_concatenated);

            % Stockage cycles bruts & outliers
            if ~isfield(CYCLES_SIGNAL, pid); CYCLES_SIGNAL.(pid)=struct(); end
            if ~isfield(CYCLES_SIGNAL.(pid), cond_ref); CYCLES_SIGNAL.(pid).(cond_ref)=struct(); end
            if ~isfield(CYCLES_SIGNAL.(pid).(cond_ref), jambes{j}); CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j})=struct(); end
            CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j}).(muscle) = all_cycles_data;

            if ~isfield(CYCLES_OUTLIERS, pid); CYCLES_OUTLIERS.(pid)=struct(); end
            if ~isfield(CYCLES_OUTLIERS.(pid), cond_ref); CYCLES_OUTLIERS.(pid).(cond_ref)=struct(); end
            if ~isfield(CYCLES_OUTLIERS.(pid).(cond_ref), jambes{j}); CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j})=struct(); end
            CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j}).(muscle) = outlier_indices;

            % Toe-off stock (une seule fois suffirait, mais on garde par muscle pour traçabilité)
            if ~isempty(all_toeoff_percentages)
                if ~isfield(CYCLES_TOEOFF, pid); CYCLES_TOEOFF.(pid)=struct(); end
                if ~isfield(CYCLES_TOEOFF.(pid), cond_ref); CYCLES_TOEOFF.(pid).(cond_ref)=struct(); end
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).mean_percentage = mean(all_toeoff_percentages);
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).std_percentage  = std(all_toeoff_percentages);
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).all_percentages = all_toeoff_percentages;
                CYCLES_TOEOFF.(pid).(cond_ref).(jambes{j}).num_cycles      = length(all_toeoff_percentages);
            end
        end

        % === Outliers globaux Plat (union tous muscles) ===
        fprintf('      Identification des outliers globaux Plat...\n');
        all_outlier_indices_plat = [];
        for m = 1:length(muscles)
            muscle = muscles{m};
            if isfield(CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j}), muscle)
                outlier_idx = CYCLES_OUTLIERS.(pid).(cond_ref).(jambes{j}).(muscle);
                all_outlier_indices_plat = union(all_outlier_indices_plat, find(outlier_idx));
            end
        end

        % === Calcul des cycles moyens BRUTS + REF_MAX_PLAT (par muscle) ===
        fprintf('      Calcul des cycles moyens BRUTS + REF_MAX_PLAT (max sur tous cycles valides)...\n');
        for m = 1:length(muscles)
            muscle = muscles{m};
            if ~isfield(CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j}), muscle), continue; end
            all_cycles = CYCLES_SIGNAL.(pid).(cond_ref).(jambes{j}).(muscle);

            % Déterminer les cycles valides pour Plat
            if ~isempty(all_outlier_indices_plat)
                total_avail = length(all_cycles);
                valid_idx   = setdiff(1:total_avail, all_outlier_indices_plat);
            else
                valid_idx   = 1:length(all_cycles);
            end
            if isempty(valid_idx), continue; end

            % Matrice [nCyclesValides x 100]
            all_points = zeros(length(valid_idx),100);
            for k = 1:length(valid_idx)
                cnum = valid_idx(k);
                if cnum <= length(all_cycles)
                    all_points(k,:) = all_cycles{cnum};
                end
            end

            % Cycle moyen brut
            cycle_moyen_brut = mean(all_points,1);

            % Référence UNIQUE : max de tous les points de tous les cycles valides
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

        % === ÉTAPE 2 — TOUTES CONDITIONS : construire SYNERGY_MATRIX + cycles normalisés avec la même ref ===
        fprintf('  === ÉTAPE 2: Construction des matrices avec normalisation unique ===\n');

        for iC = 1:length(Condition)
            cond = Condition{iC};
            fprintf('    Condition: %s\n', cond);

            % Dossier de sortie figures
            participant_folder = fullfile( ...
                'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Cycle', ...
                pid, cond);
            if ~exist(participant_folder,'dir'); mkdir(participant_folder); end

            % Si non-PLAT, extraire cycles bruts pour cette condition
            if ~strcmp(cond,'Plat')
                fprintf('      Extraction cycles + toe-off (%s)...\n', cond);
                for m = 1:length(muscles)
                    muscle = muscles{m};

                    figure_concatenated = figure('Name', sprintf('Cycles Concaténés - %s - %s - %s - %s', pid, cond, jambes{j}, muscle));
                    hold on;

                    cycle_offset   = 0;
                    cycle_positions= [];
                    muscle_maxima  = [];
                    total_cycles   = 0;
                    all_cycles_data= {};
                    all_toeoff_percentages = [];

                    for iEs = 1:length(Essai)
                        file = [pid '_' cond '_' Essai{iEs} '.c3d'];
                        try
                            data   = btkReadAcquisition(file);
                            analogs= btkGetAnalogs(data);

                            sensor_name = sensor_association.(muscle);
                            EMG_signal  = analogs.(sensor_name);

                            FreqS      = btkGetAnalogFrequency(data);
                            FreqVicon  = 100;
                            EMGfilt    = filtrage(EMG_signal, FreqS, 20, 400);
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
                                cdata = EMGenvcleaned(cycles_frames{ii});
                                cinterp = interp1(linspace(1,length(cdata),length(cdata)), cdata, linspace(1,length(cdata),100), 'pchip');

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

            % === Outliers globaux de la condition ===
            all_outlier_indices = [];
            for m = 1:length(muscles)
                muscle = muscles{m};
                if isfield(CYCLES_OUTLIERS.(pid).(cond).(jambes{j}), muscle)
                    out_idx = CYCLES_OUTLIERS.(pid).(cond).(jambes{j}).(muscle);
                    all_outlier_indices = union(all_outlier_indices, find(out_idx));
                end
            end

            % Nombre de cycles valides (prend le nb pour le 1er muscle existant)
            if isfield(CYCLES_SIGNAL.(pid).(cond).(jambes{j}), muscles{1})
                total_avail = length(CYCLES_SIGNAL.(pid).(cond).(jambes{j}).(muscles{1}));
            else
                total_avail = 0;
            end
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

            % === Construire SYNERGY_MATRIX (normalisation unique par REF_MAX_PLAT) ===
            synergy_matrix = zeros(num_valid_cycles*100, length(muscles));

            for m = 1:length(muscles)
                muscle = muscles{m};

                if ~isfield(CYCLES_SIGNAL.(pid).(cond).(jambes{j}), muscle), continue; end
                all_cycles = CYCLES_SIGNAL.(pid).(cond).(jambes{j}).(muscle);

                % Récupérer ref unique (Plat) pour ce muscle/jambe
                if isfield(REF_MAX_PLAT,pid) && isfield(REF_MAX_PLAT.(pid),jambes{j}) && isfield(REF_MAX_PLAT.(pid).(jambes{j}),muscle)
                    ref_max_plat = REF_MAX_PLAT.(pid).(jambes{j}).(muscle);
                else
                    warning('REF_MAX_PLAT manquant pour %s - %s - %s. Fallback: max individuel.', pid, jambes{j}, muscle);
                    ref_max_plat = NaN;
                end

                muscle_column = [];
                for kk = 1:length(valid_idx)
                    cnum = valid_idx(kk);
                    if cnum <= length(all_cycles)
                        cdata = all_cycles{cnum}; % 100 pts, brut (enveloppe)
                        if ~isnan(ref_max_plat) && ref_max_plat > 0
                            cdata_norm = cdata / ref_max_plat; % <<< normalisation unique
                        else
                            cmax = max(cdata); % fallback local
                            if cmax>0
                                cdata_norm = cdata / cmax;
                            else
                                cdata_norm = cdata;
                            end
                        end
                        muscle_column = [muscle_column; cdata_norm(:)];
                    end
                end

                if length(muscle_column) == num_valid_cycles*100
                    synergy_matrix(:,m) = muscle_column;
                else
                    warning('Taille incorrecte (%s - %s - %s): attendu %d, obtenu %d', muscle, pid, cond, num_valid_cycles*100, length(muscle_column));
                end

                % === Calcul cycle moyen & SD après normalisation (pour plots) ===
                n_cycles = num_valid_cycles; n_points = 100;
                if length(muscle_column) == n_cycles*n_points
                    muscle_matrix = reshape(muscle_column,[n_points,n_cycles])';
                    mean_cycle = mean(muscle_matrix,1);
                    sd_cycle   = std(muscle_matrix,0,1);

                    if ~isfield(CYCLES_MOYENS, pid), CYCLES_MOYENS.(pid)=struct(); end
                    if ~isfield(CYCLES_MOYENS.(pid), cond), CYCLES_MOYENS.(pid).(cond)=struct(); end
                    if ~isfield(CYCLES_MOYENS.(pid).(cond), jambes{j}), CYCLES_MOYENS.(pid).(cond).(jambes{j})=struct(); end
                    CYCLES_MOYENS.(pid).(cond).(jambes{j}).(muscle) = mean_cycle;

                    % Figure non affichée
                    fig = figure('Visible','off'); hold on;
                    x = linspace(0,100,n_points);
                    fill([x,fliplr(x)],[mean_cycle+sd_cycle, fliplr(mean_cycle - sd_cycle)], [1 0 0], 'FaceAlpha',0.3,'EdgeColor','none');
                    plot(x, mean_cycle, 'r','LineWidth',2);
                    title(sprintf('Cycle moyen normalisé - %s - %s - %s - %s', pid, cond, jambes{j}, muscle), 'Interpreter','none');
                    xlabel('% du cycle'); ylabel('Activation normalisée'); xlim([0 100]);

                    fig_folder = fullfile(participant_folder,'CycleMoyen');
                    if ~exist(fig_folder,'dir'); mkdir(fig_folder); end
                    fig_path = fullfile(fig_folder, sprintf('%s_%s_%s_%s_cycleMoyen.png', pid, cond, jambes{j}, muscle));
                    print(fig, fig_path, '-dpng', '-r150'); close(fig);
                end

                % === Export des cycles normalisés (même référence unique) ===
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

            end % for muscles

            % Stockage matrice de synergies
            if ~isfield(SYNERGY_MATRIX, pid), SYNERGY_MATRIX.(pid)=struct(); end
            if ~isfield(SYNERGY_MATRIX.(pid), cond), SYNERGY_MATRIX.(pid).(cond)=struct(); end
            SYNERGY_MATRIX.(pid).(cond).(jambes{j}) = synergy_matrix;

            fprintf('        SYNERGY_MATRIX: %d x %d (points concaténés x muscles)\n', size(synergy_matrix,1), size(synergy_matrix,2));
        end % for Condition
    end % for Participant
end % for jambes

% Restaurer affichage
set(0,'DefaultFigureVisible', original_visible);

% === SAUVEGARDE ===
fprintf('\n=== SAUVEGARDE DES RÉSULTATS ===\n');
save_dir  = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS';
if ~exist(save_dir,'dir'); mkdir(save_dir); end
save_file = fullfile(save_dir, [Participant{iP} '_MATRIX.mat']);

META_NORMALIZATION.method                = 'Single in-task reference (Plat)';
META_NORMALIZATION.reference_definition  = 'ref_max_plat = max over all points of all valid cycles in Plat (per participant×leg×muscle)';
META_NORMALIZATION.applied_to            = 'SYNERGY_MATRIX and CYCLES_SIGNAL_NORMALIZED for Plat/Medium/High';
META_NORMALIZATION.notes                 = 'Outliers removed before computing ref; fallback per-cycle max used only if ref missing or <=0';

save(save_file, 'SYNERGY_MATRIX', 'CYCLES_SIGNAL_NORMALIZED', ...
     'CYCLES_MOYENS_BRUTS', 'CYCLES_MOYENS', 'CYCLES_COUNT', ...
     'CYCLES_OUTLIERS', 'CYCLES_TOEOFF', 'REF_MAX_PLAT', 'META_NORMALIZATION', '-v7.3');

fprintf('Structures sauvegardées dans: %s\n', save_file);
fprintf('\n=== TRAITEMENT TERMINÉ ===\n');