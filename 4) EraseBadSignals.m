%% FILTRAGE AUTOMATIQUE DES MAUVAIS SIGNAUX (POST-NORMALISATION)
%
% OBJECTIF :
% Ce script nettoie les matrices de données en fonction d'un diagnostic 
% qualitatif consigné dans un fichier Excel ("Mapping-EMG.xlsm").
% Pour chaque muscle jugé "mauvais" (codes 1, 3 ou 80 dans Excel) :
%   1) Les données dans CYCLES_MOYENS sont remplacées par des NaN.
%   2) Les colonnes correspondantes dans SYNERGY_MATRIX sont mises à NaN.
%
% ENTRÉES :
%   - Mapping-EMG.xlsm : Fichier Excel contenant les scores de qualité (1, 2, 3, 80).
%   - Fichiers .mat originaux (issus du script 2).
%
% SORTIES :
%   - Fichiers .mat filtrés (Sauvegardés dans le dossier /FILTERED/).
%   - Log texte (*_filtrage_log.txt) résumant les suppressions par participant.
%   - SYNERGY_COLMAP : Une nouvelle structure ajoutée au .mat pour garder 
%     une trace de l'ordre des muscles en colonnes.
%
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CONFIGURATION DES CHEMINS ========================

% Dossier des résultats originaux et destination des données filtrées
matrix_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS';
output_dir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\FILTERED';
excel_path = 'C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\PhD projects\2) Projet_Surfaces_Irr\SCRIPTS\ActivationMusculaire\Mapping-EMG.xlsm';

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

% --- Paramètres de correspondance Excel ---
% Mapping des conditions avec leurs décalages de colonnes dans Excel
Condition = {'Plat',1; 'Medium',3; 'High',5}; 

% Correspondance entre le nom du muscle et sa ligne dans le fichier Excel
muscle_indices = struct( ...
    'TAprox', 13, 'TAdist', 14, 'SOL', 15, 'GM', 16, ...
    'VL', 17, 'RF', 18, 'ST', 19, 'GMED', 20);

%% ===================== LECTURE DU FICHIER QUALITÉ =======================

fprintf('Lecture du fichier Excel de mapping...\n');
try
    [~, ~, raw] = xlsread(excel_path);
    fprintf('✓ Fichier Excel lu avec succès\n');
catch ME
    error('Erreur lors de la lecture du fichier Excel: %s', ME.message);
end

% --- Identification automatique des colonnes Participants ---
participants = {}; participant_columns = {};
for col = 1:size(raw, 2)
    if ~isempty(raw{1,col}) && ischar(raw{1,col}) && startsWith(raw{1,col}, 'CTL_')
        participants{end+1} = raw{1,col}; 
        participant_columns{end+1} = col; 
    end
end

%% ===================== BOUCLE DE FILTRAGE PAR PARTICIPANT ===============

for p_idx = 1:length(participants)
    participant = participants{p_idx};
    base_col    = participant_columns{p_idx};

    % Chargement du fichier .mat original
    matrix_file = fullfile(matrix_dir, [participant '_MATRIX.mat']);
    if ~exist(matrix_file, 'file'), continue; end
    
    load(matrix_file); % Charge SYNERGY_MATRIX, CYCLES_MOYENS, etc.

    SYNERGY_COLMAP = struct(); % Pour stocker l'ordre des muscles
    removal_log_cycles = {};   % Pour le rapport final
    removal_log_synergy = {};

    for iC = 1:size(Condition,1)
        cond_name = Condition{iC,1};
        col_offset = Condition{iC,2};

        for s_idx = 1:2 % left and right
            side = if(s_idx==1, 'left', 'right');
            excel_col = base_col + col_offset + (s_idx-1);

            % --- ÉTAPE CRITIQUE : Mapping Muscles <-> Colonnes Matrix ---
            % On déduit dynamiquement quelle colonne de la matrice correspond à quel muscle
            [order, muscle2col] = derive_order_from_cycles_and_synergy( ...
                CYCLES_MOYENS, SYNERGY_MATRIX, participant, cond_name, side);

            SYNERGY_COLMAP.(participant).(cond_name).(side).order = order;
            SYNERGY_COLMAP.(participant).(cond_name).(side).muscle2col = muscle2col;

            % --- Vérification de chaque muscle via le score Excel ---
            muscle_names = fieldnames(muscle_indices);
            for m = 1:numel(muscle_names)
                muscle = muscle_names{m};
                score_excel = raw{muscle_indices.(muscle), excel_col};

                % Si le score indique un mauvais signal (1, 3 ou 80)
                if isnumeric(score_excel) && any(score_excel == [1 3 80])
                    
                    % 1. Mise à NaN dans les cycles moyens (pour les plots)
                    emg_f = ['EMG_' muscle];
                    if isfield(CYCLES_MOYENS.(participant).(cond_name), side)
                        CYCLES_MOYENS.(participant).(cond_name).(side).(emg_f)(:) = NaN;
                        removal_log_cycles{end+1} = sprintf('%s-%s-%s', cond_name, side, muscle); 
                    end

                    % 2. Mise à NaN dans la matrice de synergie (pour la NNMF)
                    if isfield(muscle2col, muscle)
                        col_idx = muscle2col.(muscle);
                        SYNERGY_MATRIX.(participant).(cond_name).(side)(:, col_idx) = NaN;
                        removal_log_synergy{end+1} = sprintf('%s-%s-col%d', cond_name, side, col_idx); 
                    end
                    fprintf('  ✗ Filtré : %s | %s | %s\n', participant, cond_name, muscle);
                end
            end
        end
    end

    % --- SAUVEGARDE DU FICHIER NETTOYÉ ---
    output_file = fullfile(output_dir, [participant '_MATRIX.mat']);
    save(output_file, 'CYCLES_COUNT','CYCLES_MOYENS','CYCLES_MOYENS_BRUTS', ...
         'CYCLES_OUTLIERS','CYCLES_SIGNAL_NORMALIZED','CYCLES_TOEOFF', ...
         'SYNERGY_MATRIX','SYNERGY_COLMAP');

    % Génération d'un log texte pour preuve de traitement
    if ~isempty(removal_log_cycles)
        write_log(output_dir, participant, removal_log_cycles, removal_log_synergy);
    end
end

fprintf('\n=== FILTRAGE TERMINÉ ===\n');

%% ===================== FONCTIONS UTILES =================================

function [order, muscle2col] = derive_order_from_cycles_and_synergy(CYCLES_MOYENS, SYNERGY_MATRIX, pid, cond, side)
    % Cette fonction assure que l'on ne supprime pas la mauvaise colonne.
    % Elle vérifie l'ordre de création des champs EMG dans la structure.
    order = {}; muscle2col = struct();
    if isfield(CYCLES_MOYENS.(pid).(cond), side)
        S = CYCLES_MOYENS.(pid).(cond).(side);
        f = fieldnames(S);
        for i = 1:numel(f)
            if strncmp(f{i},'EMG_',4), order{end+1} = erase(f{i},'EMG_'); end
        end
    end
    % Création du dictionnaire index-colonne
    for c = 1:numel(order)
        muscle2col.(order{c}) = c;
    end
end

function write_log(path, pid, log_c, log_s)
    fid = fopen(fullfile(path, [pid '_filtrage_log.txt']), 'w');
    fprintf(fid, 'LOG DE FILTRAGE - PARTICIPANT %s\n\n', pid);
    fprintf(fid, 'Muscles mis à NaN (Cycles Moyens) :\n');
    fprintf(fid, ' - %s\n', log_c{:});
    fprintf(fid, '\nColonnes mises à NaN (Synergy Matrix) :\n');
    fprintf(fid, ' - %s\n', log_s{:});
    fclose(fid);
end