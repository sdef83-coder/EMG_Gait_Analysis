%% === COHERENCE EMG-EMG (par jambe) — Script principal ===
clear; close all; clc;

addpath(genpath('C:\Users\defsil00\Documents\Script\Functions\coherence_analysis'));
cd('C:\Users\defsil00\Documents\Script\Data\jeunes_enfants'); % <- changer selon la catégorie d'âge

[num, txt, raw] = xlsread('C:\Users\defsil00\Documents\Script\Mapping-EMG.xlsm');

% === PARAMÉTRAGE UTILISATEUR ===
Participant = {['CTL_63']};   % Exemple : un seul participant

% Les numéros correspondent au décalage à droite (colonnes) dans le mapping 'raw' par rapport à 'CTL_??'
Condition = {
    'Plat',   1;
    'Medium', 3;
    'High',   5
};

Essai = {'01','02','03','04'}; % '05', '06', '07', '08', '09', '10'

raw_indices = struct( ...
    'TAprox', 13, ...
    'TAdist', 14, ...
    'SOL',    15, ...
    'GM',     16, ...
    'VL',     17, ...
    'RF',     18, ...
    'ST',     19, ...
    'GMED',   20 ...
);

alpha = 0.05; % Seuil de significativité

% === Paramètres ondelettes ===
FreqMin    = 1;
FreqMax    = 400;
Resolution = 1;
WaveNumber = 7;

FreqBands = struct( ...
    'Alpha', [8  12], ...
    'Beta',  [13  30], ...
    'Gamma', [31  60] ...
);

% Seed RNG (pour l'égalisation aléatoire des cycles)
rng(0,'twister');
seedState = rng;
fprintf('RNG seed: %d (%s)\n', seedState.Seed, seedState.Type);

% === TRAITEMENTS PAR PARTICIPANT ===
for iP = 1:numel(Participant)

    % Charge la matrice et les infos de cycles/outliers propres au participant
    load(['C:\Users\defsil00\Documents\Script\ORIGINALS\' ...
          Participant{iP} '_MATRIX.mat']);

    % Associe capteurs ↔ muscles pour left/right (défini dans ton script)
    run('C:\Users\defsil00\Documents\Script\Association.m');

    % Sécurité si la structure d’outliers est absente
    if ~exist('CYCLES_OUTLIERS','var')
        CYCLES_OUTLIERS = struct();
    end

    % Dossier de sortie
    output_dir = fullfile('C:\Users\defsil00\Documents\Script\Results\Coherence', Participant{iP});
    if ~exist(output_dir, 'dir'); mkdir(output_dir); end

    % Paires de muscles à traiter
    muscle_pairs = {
        'TAprox', 'TAdist';
        'VL',     'RF';
        'GM',     'SOL';
        'GMED',   'RF';
        'GMED',   'VL';
        'RF',     'ST'
        % ... ajouter d'autres paires si nécessaire
    };

    % Construit la table Pairs = {sensor1, sensor2, m1, m2, side, row_m1, row_m2}
    Pairs = {};
    for i = 1:size(muscle_pairs, 1)
        m1 = muscle_pairs{i,1};
        m2 = muscle_pairs{i,2};

        for side = {'left','right'}
            side_str = side{1};
            if strcmp(side_str,'left')
                assoc = sensor_association_left;
            else
                assoc = sensor_association_right;
            end

            if isfield(assoc, ['EMG_' m1]) && isfield(assoc, ['EMG_' m2]) && ...
               isfield(raw_indices, m1) && isfield(raw_indices, m2)

                Pairs(end+1,:) = { ...
                    assoc.(['EMG_' m1]), ...
                    assoc.(['EMG_' m2]), ...
                    m1, ...
                    m2, ...
                    side_str, ...
                    raw_indices.(m1), ...
                    raw_indices.(m2) ...
                };
            end
        end
    end

    % === CALCUL DES ÉVÉNEMENTS MOYENS (TO/HS) PAR CONDITION ET CÔTÉ ===
    GaitEvents = struct();

    for iC = 1:size(Condition,1)
        for main_side = {'left','right'}
            main_side_str = main_side{1};
            % Côté opposé
            if strcmp(main_side_str,'left')
                opposite_side_str = 'right';
            else
                opposite_side_str = 'left';
            end

            % Accumulateurs
            main_toe_offs_all       = [];
            opposite_toe_offs_all   = [];
            opposite_heel_strikes_all = [];

            for iEs = 1:numel(Essai)
                file = [Participant{iP} '_' Condition{iC,1} '_' Essai{iEs} '.c3d'];
                if ~isfile(file); continue; end

                data   = btkReadAcquisition(file);
                analogs = btkGetAnalogs(data);
                FreqS   = btkGetAnalogFrequency(data);
                FreqVicon = 100;

                % EMG temporaire pour les fonctions d'indexation
                flds = fieldnames(analogs);
                temp_emg = analogs.(flds{1});
                EMG_temp = [temp_emg, temp_emg];

                % Événements jambe principale
                if strcmp(main_side_str,'left')
                    main_heel_strikes = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_temp);
                    main_toe_offs     = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                else
                    main_heel_strikes = indiceRight(data, analogs, FreqS, FreqVicon, EMG_temp);
                    main_toe_offs     = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                end

                % Événements jambe opposée
                if strcmp(opposite_side_str,'left')
                    opposite_heel_strikes = indiceLeft(data, analogs, FreqS, FreqVicon, EMG_temp);
                    opposite_toe_offs     = indiceLeftTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                else
                    opposite_heel_strikes = indiceRight(data, analogs, FreqS, FreqVicon, EMG_temp);
                    opposite_toe_offs     = indiceRightTO(data, analogs, FreqS, FreqVicon, EMG_temp);
                end

                % Boucle sur les cycles
                for ii = 1:min(numel(main_heel_strikes)-1, numel(main_toe_offs))
                    if main_toe_offs(ii) > main_heel_strikes(ii) && main_toe_offs(ii) < main_heel_strikes(ii+1)
                        cycle_start  = main_heel_strikes(ii);
                        cycle_end    = main_heel_strikes(ii+1);
                        cycle_length = cycle_end - cycle_start;

                        % TO jambe principale
                        main_toeoff_percent = ((main_toe_offs(ii) - cycle_start) / cycle_length) * 100;

                        % TO jambe opposée (dans le cycle)
                        opposite_to_in_cycle = opposite_toe_offs(opposite_toe_offs > cycle_start & opposite_toe_offs < cycle_end);
                        if ~isempty(opposite_to_in_cycle)
                            opposite_toeoff_percent = ((opposite_to_in_cycle(1) - cycle_start) / cycle_length) * 100;
                        else
                            opposite_toeoff_percent = NaN;
                        end

                        % HS jambe opposée (dans le cycle)
                        opposite_hs_in_cycle = opposite_heel_strikes(opposite_heel_strikes > cycle_start & opposite_heel_strikes < cycle_end);
                        if ~isempty(opposite_hs_in_cycle)
                            opposite_heelstrike_percent = ((opposite_hs_in_cycle(1) - cycle_start) / cycle_length) * 100;
                        else
                            opposite_heelstrike_percent = NaN;
                        end

                        % Stocke si valeurs valides
                        if main_toeoff_percent > 0 && main_toeoff_percent < 100
                            main_toe_offs_all(end+1) = main_toeoff_percent; %#ok<SAGROW>
                        end
                        if ~isnan(opposite_toeoff_percent) && opposite_toeoff_percent > 0 && opposite_toeoff_percent < 100
                            opposite_toe_offs_all(end+1) = opposite_toeoff_percent; %#ok<SAGROW>
                        end
                        if ~isnan(opposite_heelstrike_percent) && opposite_heelstrike_percent > 0 && opposite_heelstrike_percent < 100
                            opposite_heel_strikes_all(end+1) = opposite_heelstrike_percent; %#ok<SAGROW>
                        end
                    end
                end
            end

            % Moyennes (ou valeurs par défaut)
            if ~isempty(main_toe_offs_all)
                GaitEvents.(Condition{iC,1}).(main_side_str).main_toeoff = mean(main_toe_offs_all);
            else
                GaitEvents.(Condition{iC,1}).(main_side_str).main_toeoff = 60;
            end

            if ~isempty(opposite_toe_offs_all)
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_toeoff = mean(opposite_toe_offs_all);
            else
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_toeoff = 10;
            end

            if ~isempty(opposite_heel_strikes_all)
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_heelstrike = mean(opposite_heel_strikes_all);
            else
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_heelstrike = 50;
            end

            fprintf('Événements moyens calculés pour %s - %s:\n', Condition{iC,1}, main_side_str);
            fprintf('  Main TO : %.1f%% | Opp TO : %.1f%% | Opp HS : %.1f%%\n', ...
                GaitEvents.(Condition{iC,1}).(main_side_str).main_toeoff, ...
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_toeoff, ...
                GaitEvents.(Condition{iC,1}).(main_side_str).opposite_heelstrike);
        end
    end

    % === TRAITEMENT PRINCIPAL (égalisation du nb total de cycles) ===
    DATA = struct();

    for iPair = 1:size(Pairs, 1)
        m1 = Pairs{iPair,3};
        m2 = Pairs{iPair,4};
        pair_label = [m1 '_' m2];

        % Trouve la colonne du participant dans 'raw'
        ii_col = [];
        for jj = 1:size(raw,2)
            if strcmp(raw{1,jj}, Participant{iP})
                ii_col = jj;
                break;
            end
        end
        if isempty(ii_col)
            warning('Participant %s introuvable dans "raw". Paire %s ignorée.', Participant{iP}, pair_label);
            continue;
        end

        % --- 1) PRE-CALCUL : indices de cycles valides/totaux par condition et côté
        ValidIdx     = struct();
        NtotalSide   = struct();
        L_total_cond = zeros(size(Condition,1),1);

        for iC = 1:size(Condition,1)
            condName = Condition{iC,1};

            for side_str_cell = {'left','right'}
                s = side_str_cell{1};
                side_offset = strcmp(s,'right'); % 0 left, 1 right

                % Vérifie validité de la paire côté s dans le mapping
                is_ok = false;
                if strcmp(s,'left')
                    assoc_side = sensor_association_left;
                else
                    assoc_side = sensor_association_right;
                end

                if isfield(assoc_side, ['EMG_' m1]) && isfield(assoc_side, ['EMG_' m2]) && ...
                   isfield(raw_indices, m1) && isfield(raw_indices, m2)
                    try
                        if raw{Pairs{iPair,6}, ii_col + Condition{iC,2} + side_offset} == 2 && ...
                           raw{Pairs{iPair,7}, ii_col + Condition{iC,2} + side_offset} == 2
                            is_ok = true;
                        end
                    catch
                        is_ok = false;
                    end
                end

                if is_ok && isfield(CYCLES_OUTLIERS, Participant{iP}) && ...
                   isfield(CYCLES_OUTLIERS.(Participant{iP}), condName) && ...
                   isfield(CYCLES_OUTLIERS.(Participant{iP}).(condName), s)

                    BadCycles = CYCLES_OUTLIERS.(Participant{iP}).(condName).(s);
                    tutu_vec  = BadCycles.(['EMG_' m1]) + BadCycles.(['EMG_' m2]); % 0 = bon sur les 2 muscles
                    ValidIdx.(condName).(s) = find(tutu_vec==0);
                    NtotalSide.(condName).(s) = numel(tutu_vec);
                else
                    ValidIdx.(condName).(s) = [];
                    NtotalSide.(condName).(s) = 0;
                end
            end

            L_total_cond(iC) = numel(ValidIdx.(condName).left) + numel(ValidIdx.(condName).right);
        end

        % --- 2) CIBLE = min(total gauche+droite) parmi les conditions (strictement positif)
pos_mask = L_total_cond > 0;
if any(pos_mask)
    L_target_eff = min(L_total_cond(pos_mask));   % minimum strictement positif
else
    L_target_eff = 0;  % aucune condition n'a de cycles valides
end

% message d'info uniquement si TOUT est nul
if L_target_eff == 0
    fprintf('⚠ Paire %s : aucune condition n''a de cycles valides (%s). Paire ignorée.\n', pair_label, Participant{iP});
end

% Stocker une cible dans meta (par condition) pour traçabilité
for iiC = 1:size(Condition,1)
    condName_ii = Condition{iiC,1};
    DATA.(condName_ii).meta.(['L_target_' m1 '_' m2]) = L_target_eff;
    if isfield(NtotalSide, condName_ii)
        DATA.(condName_ii).meta.(['Ltot_' m1 '_' m2 '_left'])  = NtotalSide.(condName_ii).left;
        DATA.(condName_ii).meta.(['Ltot_' m1 '_' m2 '_right']) = NtotalSide.(condName_ii).right;
    end
end

% Si tout est nul, on passera à la suite mais rien ne sera retenu

     % --- 3) Construit les masques de cycles à garder (KeepMask)
KeepMask = struct();
for iC = 1:size(Condition,1)
    condName = Condition{iC,1};
    nL = NtotalSide.(condName).left;
    nR = NtotalSide.(condName).right;

    KeepMask.(condName).left  = false(nL,1);
    KeepMask.(condName).right = false(nR,1);

    idxL = ValidIdx.(condName).left;
    idxR = ValidIdx.(condName).right;

    union_sides = [
        [repmat({'left'},  numel(idxL), 1), num2cell(idxL(:))];
        [repmat({'right'}, numel(idxR), 1), num2cell(idxR(:))]
    ];
    n_union = size(union_sides,1);

    if n_union == 0
        % Rien à garder pour cette condition
        continue;
    end

    if L_target_eff == 0
        % Aucune condition globale non nulle -> garder tout ce qui est valide (pas d'égalisation)
        if ~isempty(idxL), KeepMask.(condName).left(idxL)  = true; end
        if ~isempty(idxR), KeepMask.(condName).right(idxR) = true; end
    else
        % Égalisation sur le minimum strictement positif
        if n_union <= L_target_eff
            if ~isempty(idxL), KeepMask.(condName).left(idxL)  = true; end
            if ~isempty(idxR), KeepMask.(condName).right(idxR) = true; end
        else
            sel = randperm(n_union, L_target_eff);
            chosen = union_sides(sel,:);
            for k = 1:size(chosen,1)
                if strcmp(chosen{k,1}, 'left')
                    KeepMask.(condName).left(chosen{k,2}) = true;
                else
                    KeepMask.(condName).right(chosen{k,2}) = true;
                end
            end
        end
    end
end

        % --- 4) Calculs spectraux / cohérence en ne gardant que KeepMask
        for iC = 1:size(Condition,1)
            condName = Condition{iC,1};

            % Côté porté par la ligne Pairs (pas de boucle left/right ici)
            side_str = Pairs{iPair,5};
            Side = {side_str, double(strcmp(side_str,'right'))}; % {'left',0} ou {'right',1}

            % Validité de la paire (colonnes mapping)
            is_ok = false;
            try
                if raw{Pairs{iPair,6}, ii_col + Condition{iC,2} + Side{1,2}} == 2 && ...
                   raw{Pairs{iPair,7}, ii_col + Condition{iC,2} + Side{1,2}} == 2
                    is_ok = true;
                end
            catch
                is_ok = false;
            end
            if ~is_ok, continue; end

            % Indices "mauvais cycles" pour CE côté
            if isfield(CYCLES_OUTLIERS, Participant{iP}) && ...
               isfield(CYCLES_OUTLIERS.(Participant{iP}), condName) && ...
               isfield(CYCLES_OUTLIERS.(Participant{iP}).(condName), side_str)
                BadCycles = CYCLES_OUTLIERS.(Participant{iP}).(condName).(side_str);
                tutu_vec  = BadCycles.(['EMG_' m1]) + BadCycles.(['EMG_' m2]); % 0 = à garder
            else
                tutu_vec  = [];
            end

            % Accumulateurs (dimensionnés dynamiquement au 1er cycle retenu)
            PowSpec_s1    = [];
            PowSpec_s2    = [];
            cross_spectrum = [];
            NcycleP = 0; % index de cycle "global" pour faire correspondre tutor/KeepMask
            NcycleC = 0; % nombre de cycles effectivement conservés

            for iEs = 1:numel(Essai)
                file = [Participant{iP} '_' condName '_' Essai{iEs} '.c3d'];
                if ~isfile(file), continue; end

                data    = btkReadAcquisition(file);
                analogs = btkGetAnalogs(data);

                % EMG de la paire
                EMG = [];
                EMG(:,1) = analogs.(Pairs{iPair,1});
                EMG(:,2) = analogs.(Pairs{iPair,2});

                FreqS     = btkGetAnalogFrequency(data);
                FreqVicon = 100; %#ok<NASGU> % (non utilisé ici directement)

                % Cycles du bon côté
                if strcmp(Side{1},'left')
                    Cycles = indiceLeft(data, analogs, FreqS, FreqVicon, EMG);
                else
                    Cycles = indiceRight(data, analogs, FreqS, FreqVicon, EMG);
                end

                % Ondelette + prétraitements
                Args = WaveletParameters(FreqMin, FreqMax, Resolution, WaveNumber, FreqS);
                EMG  = filtrage(EMG, FreqS, 8, 400);
                EMG  = EMG - mean(EMG,1);

                [TFR(:,:,1), period, ~, ~] = wavelet(EMG(:,1), Args.DT, Args.Pad, Args.DJ, Args.S0, Args.J1, Args.Mother, Args.Cycles);
                [TFR(:,:,2), ~,     ~, ~]  = wavelet(EMG(:,2), Args.DT, Args.Pad, Args.DJ, Args.S0, Args.J1, Args.Mother, Args.Cycles);
                Freq = 1 ./ period;

                % Parcours des cycles de l'essai
                for iCycles = 1:(numel(Cycles)-1)
                    NcycleP = NcycleP + 1; % index global

                    % 1) cycle valide (tutu==0) OU pas de vector tutor ; 2) cycle sélectionné par KeepMask
                    is_valid = (NcycleP > numel(tutu_vec)) || (tutu_vec(NcycleP)==0);
                    in_mask  = (NcycleP <= numel(KeepMask.(condName).(side_str))) && KeepMask.(condName).(side_str)(NcycleP);

                    if is_valid && in_mask
                        % Découpe + interpolation 1000 points
                        TFR_cycle(:,:,1) = TFR(:, Cycles(iCycles):Cycles(iCycles+1), 1);
                        TFR_cycle(:,:,2) = TFR(:, Cycles(iCycles):Cycles(iCycles+1), 2);

                        [X, Y]  = meshgrid(1:size(TFR_cycle,2), 1:size(TFR_cycle,1));
                        [Xq, Yq]= meshgrid(linspace(1,size(TFR_cycle,2),1000), 1:size(TFR_cycle,1));

                        TFR_int(:,:,1) = interp2(X, Y, TFR_cycle(:,:,1), Xq, Yq, 'spline');
                        TFR_int(:,:,2) = interp2(X, Y, TFR_cycle(:,:,2), Xq, Yq, 'spline');

                        % Initialise les accumulateurs à la bonne taille au 1er cycle retenu
                        if isempty(PowSpec_s1)
                            PowSpec_s1     = zeros(size(TFR_int,1), size(TFR_int,2));
                            PowSpec_s2     = zeros(size(TFR_int,1), size(TFR_int,2));
                            cross_spectrum = zeros(size(TFR_int,1), size(TFR_int,2));
                        end

                        % Accumulations (sommes)
                        PowSpec_s1     = PowSpec_s1     + abs(TFR_int(:,:,1)).^2;
                        PowSpec_s2     = PowSpec_s2     + abs(TFR_int(:,:,2)).^2;
                        cross_spectrum = cross_spectrum + (TFR_int(:,:,1)) .* conj(TFR_int(:,:,2));

                        NcycleC = NcycleC + 1;

                        clear TFR_int TFR_cycle
                    end
                end

                clear TFR
            end % essais

            if NcycleC == 0
                warning('Aucun cycle retenu après égalisation | %s-%s | %s | %s', m1, m2, condName, side_str);
                continue;
            end

            % === Cohérence ===
            Pxx_sum = PowSpec_s1;
            Pyy_sum = PowSpec_s2;
            Pxy_sum = cross_spectrum;
            L_cycles = NcycleC;

            PowSpec_s1_mean    = PowSpec_s1 / NcycleC;
            PowSpec_s2_mean    = PowSpec_s2 / NcycleC;
            cross_spectrum_mean = cross_spectrum / NcycleC;

            Coherence = abs(cross_spectrum_mean).^2 ./ (PowSpec_s1_mean .* PowSpec_s2_mean);

            % === Événements & sous-phases ===
            ge = GaitEvents.(condName).(side_str);
            main_toeoff_index      = max(1,   min(round(ge.main_toeoff     * 10), 1000));
            opposite_toeoff_index  = max(1,   min(round(ge.opposite_toeoff * 10), max(1,main_toeoff_index-1)));
            opposite_heelstrike_index = max(opposite_toeoff_index+1, min(round(ge.opposite_heelstrike*10), max(2,main_toeoff_index-1)));

            loading_response_indices = 1:opposite_toeoff_index;
            midstance_indices        = (opposite_toeoff_index+1):opposite_heelstrike_index;
            preswing_indices         = (opposite_heelstrike_index+1):main_toeoff_index;
            swing_phase_indices      = (main_toeoff_index+1):1000;

            % === Moyennes par bande (Alpha/Beta/Gamma) ===
            mean_coherence_loading  = struct(); 
            mean_coherence_midstance= struct();
            mean_coherence_preswing = struct();
            mean_coherence_swing    = struct();
            mean_coherence_full     = struct();

            for bandName = fieldnames(FreqBands)'
                band  = bandName{1};
                range = FreqBands.(band);
                idx_band = find(Freq >= range(1) & Freq <= range(2));

                mean_coherence_full.(band) = mean(mean(Coherence(idx_band,:), 2), 'omitnan');

                if ~isempty(loading_response_indices)
                    mean_coherence_loading.(band)   = mean(mean(Coherence(idx_band, loading_response_indices), 2), 'omitnan');
                else
                    mean_coherence_loading.(band)   = NaN;
                end

                if ~isempty(midstance_indices)
                    mean_coherence_midstance.(band) = mean(mean(Coherence(idx_band, midstance_indices), 2), 'omitnan');
                else
                    mean_coherence_midstance.(band) = NaN;
                end

                if ~isempty(preswing_indices)
                    mean_coherence_preswing.(band)  = mean(mean(Coherence(idx_band, preswing_indices), 2), 'omitnan');
                else
                    mean_coherence_preswing.(band)  = NaN;
                end

                if ~isempty(swing_phase_indices)
                    mean_coherence_swing.(band)     = mean(mean(Coherence(idx_band, swing_phase_indices), 2), 'omitnan');
                else
                    mean_coherence_swing.(band)     = NaN;
                end

                % Enregistrement dans DATA
                DATA.(condName).(side_str).(['MeanCoherence_' band '_' m1 '_' m2])                    = mean_coherence_full.(band);
                DATA.(condName).(side_str).(['MeanCoherence_LoadingResponse_' band '_' m1 '_' m2])    = mean_coherence_loading.(band);
                DATA.(condName).(side_str).(['MeanCoherence_MidStance_' band '_' m1 '_' m2])          = mean_coherence_midstance.(band);
                DATA.(condName).(side_str).(['MeanCoherence_PreSwing_' band '_' m1 '_' m2])           = mean_coherence_preswing.(band);
                DATA.(condName).(side_str).(['MeanCoherence_Swing_' band '_' m1 '_' m2])              = mean_coherence_swing.(band);
            end

            % Seuil de significativité
            if NcycleC <= 1
                seuil = NaN;
            else
                seuil = 1 - alpha^(1/(NcycleC - 1));
            end

            % === Sauvegardes ===
            DATA.(condName).(side_str).(['GaitEvents_' m1 '_' m2]) = struct( ...
                'main_toeoff',        ge.main_toeoff, ...
                'opposite_toeoff',    ge.opposite_toeoff, ...
                'opposite_heelstrike',ge.opposite_heelstrike);

            DATA.(condName).(side_str).(['Coherence_LoadingResponse_' m1 '_' m2]) = Coherence(:, loading_response_indices);
            DATA.(condName).(side_str).(['Coherence_MidStance_'       m1 '_' m2]) = Coherence(:, midstance_indices);
            DATA.(condName).(side_str).(['Coherence_PreSwing_'        m1 '_' m2]) = Coherence(:, preswing_indices);
            DATA.(condName).(side_str).(['Coherence_Swing_'           m1 '_' m2]) = Coherence(:, swing_phase_indices);

            DATA.(condName).(side_str).(['LoadingResponseIndices_' m1 '_' m2]) = loading_response_indices;
            DATA.(condName).(side_str).(['MidStanceIndices_'       m1 '_' m2]) = midstance_indices;
            DATA.(condName).(side_str).(['PreSwingIndices_'        m1 '_' m2]) = preswing_indices;
            DATA.(condName).(side_str).(['SwingIndices_'           m1 '_' m2]) = swing_phase_indices;

            DATA.(condName).(side_str).(['PowSpec' m1]) = PowSpec_s1;
            DATA.(condName).(side_str).(['PowSpec' m2]) = PowSpec_s2;

            % Pour compatibilité : cross-spectrum moyen au carré (module^2)
            DATA.(condName).(side_str).(['Cross_Spectrum_' m1 '_' m2]) = abs(cross_spectrum / NcycleC).^2;

            DATA.(condName).(side_str).(['Coherence_' m1 '_' m2]) = Coherence;
            DATA.(condName).(side_str).(['Ncycle_'    m1 '_' m2]) = NcycleC;

            % Métadonnées
            DATA.(condName).meta.(['L_target_eff_'           m1 '_' m2]) = L_target_eff;
            DATA.(condName).meta.(['Ltot_'               m1 '_' m2 '_left'])  = NtotalSide.(condName).left;
            DATA.(condName).meta.(['Ltot_'               m1 '_' m2 '_right']) = NtotalSide.(condName).right;
            DATA.(condName).meta.RNG.Type  = seedState.Type;
            DATA.(condName).meta.RNG.Seed  = seedState.Seed;

            DATA.(condName).(side_str).(['Seuil_' m1 '_' m2]) = seuil;
            DATA.(condName).(side_str).(['PxxSum_' m1 '_' m2]) = Pxx_sum;
            DATA.(condName).(side_str).(['PyySum_' m1 '_' m2]) = Pyy_sum;
            DATA.(condName).(side_str).(['PxySum_' m1 '_' m2]) = Pxy_sum; % complexe
            DATA.(condName).(side_str).(['L_' m1 '_' m2])      = L_cycles;
            DATA.(condName).(side_str).Freq = Freq;

            % === Figures ===
            Time = linspace(0,100,1000);
            muscle1_name = m1;
            muscle2_name = m2;

            figure('Position',[100,100,1400,800]);

            subplot(2,4,1);
            imagesc(Time, Freq, PowSpec_s1);
            title(['Power Spec ' muscle1_name]); xlabel('% Cycle'); ylabel('Fréquence (Hz)'); colorbar;

            subplot(2,4,2);
            imagesc(Time, Freq, PowSpec_s2);
            title(['Power Spec ' muscle2_name]); xlabel('% Cycle'); ylabel('Fréquence (Hz)'); colorbar;

            subplot(2,4,3);
            imagesc(Time, Freq, Coherence > seuil);
            title([muscle1_name '-' muscle2_name ' Coherence']); hold on;
            line([ge.opposite_toeoff ge.opposite_toeoff], [min(Freq) max(Freq)], 'Color','blue',  'LineWidth',2);
            line([ge.opposite_heelstrike ge.opposite_heelstrike], [min(Freq) max(Freq)], 'Color','green', 'LineWidth',2);
            line([ge.main_toeoff ge.main_toeoff], [min(Freq) max(Freq)], 'Color','red',   'LineWidth',2);
            xlabel('% Cycle'); ylabel('Fréquence (Hz)');

            subplot(2,4,4);
            plot(Freq, mean(Coherence,2)); hold on;
            line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Cycle complet'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            subplot(2,4,5);
            if ~isempty(loading_response_indices), plot(Freq, mean(Coherence(:,loading_response_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Loading Response'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            subplot(2,4,6);
            if ~isempty(midstance_indices), plot(Freq, mean(Coherence(:,midstance_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Mid Stance'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            subplot(2,4,7);
            if ~isempty(preswing_indices), plot(Freq, mean(Coherence(:,preswing_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Pre-swing'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            subplot(2,4,8);
            if ~isempty(swing_phase_indices), plot(Freq, mean(Coherence(:,swing_phase_indices),2)); end
            hold on; line([0 400],[seuil seuil],'Color','red');
            xlim([1 400]); ylim([0 1]);
            title('Cohérence - Swing'); xlabel('Fréquence (Hz)'); ylabel('Cohérence');

            sgtitle(['Condition : ' condName ' | Côté : ' side_str ' | ' muscle1_name ' - ' muscle2_name]);

            fig_name = sprintf('%s_%s_%s-%s', condName, side_str, muscle1_name, muscle2_name);
            saveas(gcf, fullfile(output_dir, [fig_name '.png']));
            close(gcf);

            fprintf('\n=== Résultats pour %s-%s | %s | %s ===\n', muscle1_name, muscle2_name, condName, side_str);
            fprintf('  Opp TO : %.1f%% | Opp HS : %.1f%% | Main TO : %.1f%%\n', ...
                ge.opposite_toeoff, ge.opposite_heelstrike, ge.main_toeoff);

        end % Conditions
    end % Pairs

    % Sauvegarde finale
    save(fullfile(output_dir, ['Coherence_' Participant{iP} '.mat']), 'DATA');
    fprintf('Script terminé avec succès pour %s !\n', Participant{iP});

end % Participants