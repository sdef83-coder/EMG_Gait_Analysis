%% FILTRATION AUTOMATIQUE DES MAUVAIS SIGNAUX (CYCLES_MOYENS + SYNERGY_MATRIX)
% Déduction dynamique de l'ordre des muscles depuis CYCLES_MOYENS (champs EMG_*).
% Sauvegarde un mapping SYNERGY_COLMAP (order + muscle2col) et un log texte.

clc; clear; close all;

% --- CHEMINS (à adapter si besoin) ---
excel_path = 'C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\PhD projects\2) Projet_Surfaces_Irr\SCRIPTS\ActivationMusculaire\Mapping-EMG.xlsm';
matrix_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS';
output_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\FILTERED';

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% --- Conditions (comme dans ton code) ---
Condition = {'Plat',1; 'Medium',3; 'High',5};  % décalages Excel

% --- Lignes Excel associées aux muscles dans Mapping-EMG.xlsm ---
muscle_indices = struct( ...
    'TAprox', 13, ...
    'TAdist', 14, ...
    'SOL',   15, ...
    'GM',    16, ...
    'VL',    17, ...
    'RF',    18, ...
    'ST',    19, ...
    'GMED',  20);

% --- Lecture Excel ---
fprintf('Lecture du fichier Excel de mapping...\n');
try
    [num, txt, raw] = xlsread(excel_path); %#ok<ASGLU>
    fprintf('✓ Fichier Excel lu avec succès\n');
catch ME
    error('Erreur lors de la lecture du fichier Excel: %s', ME.message);
end

% --- Identifier les participants (entêtes "CTL_..") ---
participants = {};
participant_columns = {};
for col = 1:size(raw, 2)
    if ~isempty(raw{1,col}) && ischar(raw{1,col})
        cell_content = raw{1,col};
        if startsWith(cell_content, 'CTL_')
            participants{end+1} = cell_content; %#ok<SAGROW>
            participant_columns{end+1} = col;   %#ok<SAGROW>
        end
    end
end
if isempty(participants)
    warning('Aucun participant "CTL_*" trouvé dans la 1ère ligne de %s.', excel_path);
else
    fprintf('Participants trouvés: %s\n', strjoin(participants, ', '));
end

% =========================================================================
% == TRAITEMENT PAR PARTICIPANT ===========================================
% =========================================================================
for p_idx = 1:length(participants)
    participant = participants{p_idx};
    base_col   = participant_columns{p_idx};

    fprintf('\n=== Traitement de %s ===\n', participant);

    % Charger la matrice du participant
    matrix_file = fullfile(matrix_dir, [participant '_MATRIX.mat']);
    if ~exist(matrix_file, 'file')
        warning('Fichier matrice non trouvé pour %s: %s', participant, matrix_file);
        continue;
    end
    fprintf('Chargement de la matrice: %s\n', matrix_file);
    load(matrix_file, 'CYCLES_COUNT','CYCLES_MOYENS','CYCLES_MOYENS_BRUTS','CYCLES_OUTLIERS','CYCLES_SIGNAL_NORMALIZED','CYCLES_TOEOFF','SYNERGY_MATRIX');

    % Initialiser la structure de mapping des colonnes SYNERGY
    SYNERGY_COLMAP = struct();

    total_removed_cycles = 0;  % nb de signaux mis à NaN dans CYCLES_MOYENS
    total_nan_synergy    = 0;  % nb de colonnes mises à NaN dans SYNERGY_MATRIX
    removal_log_cycles   = {};
    removal_log_synergy  = {};

    sides = {'left','right'};
    side_offsets = [0,1]; % left=col+0 ; right=col+1

    % Parcourir conditions
    for iC = 1:size(Condition,1)
        condition_name = Condition{iC,1};
        col_offset     = Condition{iC,2};

        for s_idx = 1:numel(sides)
            side = sides{s_idx};
            side_offset = side_offsets(s_idx);

            % Colonne Excel à lire pour ce participant/condition/side
            excel_col = base_col + col_offset + side_offset;

            % --- Déduire dynamiquement l'ordre des muscles pour SYNERGY depuis CYCLES_MOYENS ---
            [order, muscle2col] = derive_order_from_cycles_and_synergy( ...
                CYCLES_MOYENS, SYNERGY_MATRIX, participant, condition_name, side);

            % Ranger cette correspondance
            if ~isfield(SYNERGY_COLMAP, participant); SYNERGY_COLMAP.(participant) = struct(); end
            if ~isfield(SYNERGY_COLMAP.(participant), condition_name); SYNERGY_COLMAP.(participant).(condition_name) = struct(); end
            SYNERGY_COLMAP.(participant).(condition_name).(side).order = order;
            SYNERGY_COLMAP.(participant).(condition_name).(side).muscle2col = muscle2col;

            % --- Boucle muscles (pilotée par muscle_indices pour lire Excel) ---
            muscle_names = fieldnames(muscle_indices);
            for m_idx = 1:numel(muscle_names)
                muscle = muscle_names{m_idx};
                muscle_row = muscle_indices.(muscle);
                emg_field = ['EMG_' muscle];

                if excel_col <= size(raw,2) && muscle_row <= size(raw,1)
                    excel_value = raw{muscle_row, excel_col};

                    % -------- Filtrage CYCLES_MOYENS (ton code existant) --------
                    if isnumeric(excel_value) && any(excel_value == [1 3 80]) % mauvais signal
                        if isfield(CYCLES_MOYENS, participant) && ...
                           isfield(CYCLES_MOYENS.(participant), condition_name) && ...
                           isfield(CYCLES_MOYENS.(participant).(condition_name), side) && ...
                           isfield(CYCLES_MOYENS.(participant).(condition_name).(side), emg_field)

                            CYCLES_MOYENS.(participant).(condition_name).(side).(emg_field)(:) = NaN;
                            total_removed_cycles = total_removed_cycles + 1;
                            removal_log_cycles{end+1} = sprintf('%s - %s - %s - %s', participant, condition_name, side, emg_field);
                            fprintf('  ✗ CYCLES_MOYENS filtré: %s/%s/%s/%s (Excel=%s)\n', participant, condition_name, side, emg_field, mat2str(excel_value));
                        else
                            fprintf('  ⚠ Structure CYCLES_MOYENS absente: %s/%s/%s/%s\n', participant, condition_name, side, emg_field);
                        end

                        % -------- Filtrage SYNERGY_MATRIX (colonne -> NaN) -----------
                        if isfield(SYNERGY_MATRIX, participant) && ...
                           isfield(SYNERGY_MATRIX.(participant), condition_name) && ...
                           isfield(SYNERGY_MATRIX.(participant).(condition_name), side)

                            M = SYNERGY_MATRIX.(participant).(condition_name).(side);
                            if isnumeric(M) && ~isempty(M)
                                if isfield(muscle2col, muscle)
                                    col_idx = muscle2col.(muscle);
                                    if col_idx <= size(M,2)
                                        M(:, col_idx) = NaN;
                                        SYNERGY_MATRIX.(participant).(condition_name).(side) = M; % ré-injecter
                                        total_nan_synergy = total_nan_synergy + 1;
                                        removal_log_synergy{end+1} = sprintf('%s - %s - %s - col %d (%s)', participant, condition_name, side, col_idx, muscle);
                                        fprintf('  ✗ SYNERGY_MATRIX NaN colonne %d (%s/%s/%s)\n', col_idx, participant, condition_name, side);
                                    else
                                        fprintf('  ⚠ Index colonne > nb colonnes (muscle=%s) dans SYNERGY: %s/%s/%s\n', muscle, participant, condition_name, side);
                                    end
                                else
                                    fprintf('  ⚠ Muscle "%s" absent de l ordre détecté (%s/%s/%s)\n', muscle, participant, condition_name, side);
                                end
                            else
                                fprintf('  ⚠ SYNERGY_MATRIX vide ou non num. (%s/%s/%s)\n', participant, condition_name, side);
                            end
                        else
                            fprintf('  ⚠ Chemin SYNERGY_MATRIX inexistant (%s/%s/%s)\n', participant, condition_name, side);
                        end

                    elseif isnumeric(excel_value) && excel_value == 2
                        fprintf('  ✓ Conservé: %s/%s/%s/%s (Excel=%s)\n', participant, condition_name, side, emg_field, mat2str(excel_value));
                    else
                        fprintf('  ? Valeur inconnue Excel: %s/%s/%s/%s (Excel=%s)\n', participant, condition_name, side, emg_field, mat2str(excel_value));
                    end
                end
            end
        end
    end

    % --- Sauvegarde (avec SYNERGY_COLMAP ajouté) ---
    output_file = fullfile(output_dir, [participant '_MATRIX.mat']);
    save(output_file, 'CYCLES_COUNT','CYCLES_MOYENS','CYCLES_MOYENS_BRUTS','CYCLES_OUTLIERS','CYCLES_SIGNAL_NORMALIZED','CYCLES_TOEOFF','SYNERGY_MATRIX','SYNERGY_COLMAP');
    fprintf('\n✓ Matrice sauvegardée (FILTERED): %s\n', output_file);
    fprintf('Total filtré CYCLES_MOYENS: %d | Colonnes NaN SYNERGY: %d\n', total_removed_cycles, total_nan_synergy);

    % --- Logs texte optionnels ---
    if ~isempty(removal_log_cycles) || ~isempty(removal_log_synergy)
        log_file = fullfile(output_dir, [participant '_filtrage_log.txt']);
        fid = fopen(log_file, 'w');
        fprintf(fid, 'Filtrage pour %s\n\n', participant);
        if ~isempty(removal_log_cycles)
            fprintf(fid, 'CYCLES_MOYENS (mis à NaN):\n');
            for i=1:numel(removal_log_cycles), fprintf(fid, ' - %s\n', removal_log_cycles{i}); end
            fprintf(fid, '\n');
        end
        if ~isempty(removal_log_synergy)
            fprintf(fid, 'SYNERGY_MATRIX (colonnes NaN):\n');
            for i=1:numel(removal_log_synergy), fprintf(fid, ' - %s\n', removal_log_synergy{i}); end
        end
        fclose(fid);
        fprintf('Log de filtrage: %s\n', log_file);
    end
end

fprintf('\n=== Filtration terminée (CYCLES_MOYENS + SYNERGY_MATRIX) ===\n');

% -------------------------------------------------------------------------
% UTILITAIRE: déduit l'ordre des colonnes via CYCLES_MOYENS (fallback: champs EMG_)
% -------------------------------------------------------------------------
function [order, muscle2col] = derive_order_from_cycles_and_synergy(CYCLES_MOYENS, SYNERGY_MATRIX, participant, condition_name, side)
% Renvoie:
%  - order: cellstr ordre des muscles détectés
%  - muscle2col: struct muscle -> index de colonne
%
% Stratégie:
% 1) Si CYCLES_MOYENS.(...).(side) contient 'MUSCLE_ORDER' (cellstr) ou 'order',
%    on l'utilise en priorité.
% 2) Sinon, on énumère les champs 'EMG_*' dans l'ordre retourné par fieldnames
%    (>= R2018b: ordre de création, ce qui matche généralement ton pipeline).
% 3) On tronque/ajuste selon le nb de colonnes réellement présent dans SYNERGY_MATRIX.

order = {};
muscle2col = struct();

% 1) Candidates depuis CYCLES_MOYENS
ord_candidates = {};
if isfield(CYCLES_MOYENS, participant) && ...
   isfield(CYCLES_MOYENS.(participant), condition_name) && ...
   isfield(CYCLES_MOYENS.(participant).(condition_name), side)

    S = CYCLES_MOYENS.(participant).(condition_name).(side);

    if isfield(S, 'MUSCLE_ORDER') && iscellstr(S.MUSCLE_ORDER)
        ord_candidates = S.MUSCLE_ORDER(:).';
    elseif isfield(S, 'order') && iscellstr(S.order)
        ord_candidates = S.order(:).';
    else
        f = fieldnames(S);
        for i = 1:numel(f)
            fn = f{i};
            if strncmp(fn,'EMG_',4)
                ord_candidates{end+1} = erase(fn,'EMG_'); %#ok<AGROW>
            end
        end
    end
end

% 2) Nb de colonnes dans SYNERGY_MATRIX
ncol = NaN;
if isfield(SYNERGY_MATRIX, participant) && ...
   isfield(SYNERGY_MATRIX.(participant), condition_name) && ...
   isfield(SYNERGY_MATRIX.(participant).(condition_name), side)
    M = SYNERGY_MATRIX.(participant).(condition_name).(side);
    if isnumeric(M) && ~isempty(M)
        ncol = size(M,2);
    end
end

% 3) Construire 'order'
if isempty(ord_candidates)
    order = {};
    return;
end

if ~isnan(ncol)
    order = ord_candidates(1:min(numel(ord_candidates), ncol));
else
    order = ord_candidates; % si SYNERGY absent/vide, on retourne la liste issue de CYCLES_MOYENS
end

% 4) Mapping muscle -> colonne
for c = 1:numel(order)
    m = order{c};
    % éviter noms vides/invalides
    if ischar(m) && ~isempty(m)
        m = strrep(m, ' ', ''); % hygiène minimale
        muscle2col.(m) = c;
    end
end
end