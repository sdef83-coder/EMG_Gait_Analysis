%% REGROUPEMENT DES DONNÉES (MERGE) POUR ANALYSES STATISTIQUES ET SYNERGIES
%
% OBJECTIF :
% Ce script fusionne les données de tous les participants d'un groupe (ex: Enfants, Adultes).
% Il génère trois types de fichiers de sortie :
%   1) SPM-EMG-XXX.mat : Pour les analyses statistiques de type Statistical Parametric Mapping.
%   2) SYNERGY-XXX.mat : Matrices globales concaténées pour l'analyse NNMF.
%   3) PARTICIPANT-CYCLES-XXX.mat : Suivi individuel des cycles par participant.
%
% ENTRÉES :
%   - Fichiers .mat individuels (provenant du dossier /FILTERED/ ou /ORIGINALS/).
%
% SORTIES :
%   - Structures fusionnées par condition (Plat, Medium, High).
%   - Affichage console des valeurs max tous les 10% du cycle (contrôle qualité).
%
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CONFIGURATION DU GROUPE D'ÉTUDE ==================

% Dossier racine des résultats
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix')

% Paramètres du groupe (À MODIFIER SELON L'ANALYSE)
groupe_a_etudier = 'jeunes_enfants'; % Utilisé pour le nom du fichier de sortie
Participant = {'CTL_14', 'CTL_15', 'CTL_23', 'CTL_37', 'CTL_40', 'CTL_53', 'CTL_63'}; 

% Paramètres fixes
Conditions = {'Plat', 'Medium', 'High'};  
Muscles    = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'};
Jambes     = {'left', 'right'};

%% ===================== INITIALISATION DES STRUCTURES ====================

% Stockage pour SPM (Cycles moyens)
mean_cycles_combined = struct(); 
% Stockage pour Synergies (Matrices concaténées)
SYNERGY_Adulte = struct(); 
% Suivi des cycles
normalized_participant_cycles = struct();
participant_tracking = struct();

fprintf('Fusion des données pour le groupe : %s\n', groupe_a_etudier);

%% ===================== BOUCLE DE FUSION DES DONNÉES =====================

for iC = 1:length(Conditions)
    condition = Conditions{iC};
    fprintf('\n--- Traitement Condition : %s ---\n', condition);
    
    condition_matrix = []; % Matrice pour concaténer tous les cycles du groupe
    muscle_names     = {}; 

    for iP = 1:length(Participant)
        pid = Participant{iP};
        
        % Chargement du fichier filtré du participant
        filename = fullfile('FILTERED', [pid '_MATRIX.mat']);
        if ~exist(filename, 'file')
            warning('Fichier manquant pour %s : %s', pid, filename);
            continue;
        end
        load(filename); % Charge SYNERGY_MATRIX, CYCLES_MOYENS, etc.

        for iJ = 1:length(Jambes)
            leg = Jambes{iJ};
            
            % --- 1. Extraction pour SPM (Cycles Moyens) ---
            if isfield(CYCLES_MOYENS.(pid), condition) && isfield(CYCLES_MOYENS.(pid).(condition), leg)
                for iM = 1:length(Muscles)
                    m_name = Muscles{iM};
                    if isfield(CYCLES_MOYENS.(pid).(condition).(leg), m_name)
                        % Stockage du cycle moyen (1x100 points)
                        mean_cycles_combined.(condition).(leg).(m_name)(iP, :) = ...
                            CYCLES_MOYENS.(pid).(condition).(leg).(m_name);
                    end
                end
            end

            % --- 2. Extraction pour Synergies (Matrices de points) ---
            % On récupère la matrice de synergie (N cycles * 100 pts x M muscles)
            if isfield(SYNERGY_MATRIX.(pid), condition) && isfield(SYNERGY_MATRIX.(pid).(condition), leg)
                current_matrix = SYNERGY_MATRIX.(pid).(condition).(leg);
                
                % Concaténation verticale (tous les participants les uns sous les autres)
                condition_matrix = [condition_matrix; current_matrix]; %#ok<AGROW>
                
                % Stockage des métadonnées (nom des muscles) la première fois
                if isempty(muscle_names) && isfield(SYNERGY_COLMAP.(pid).(condition), leg)
                    muscle_names = SYNERGY_COLMAP.(pid).(condition).(leg).order;
                end
            end
        end
    end

    % ===================== ANALYSE DES VALEURS PAR 10% ==================
    % Ce bloc permet de vérifier rapidement la distribution de l'activation
    % sur le cycle de marche (tous les 10%) pour le groupe.
    
    fprintf('Récapitulatif de l''activation max (tous les 10%%) :\n');
    for iM = 1:length(Muscles)
        m_name = Muscles{iM};
        if isfield(mean_cycles_combined.(condition).left, m_name)
            data_avg = mean(mean_cycles_combined.(condition).left.(m_name), 1, 'omitnan');
            % Extraction des valeurs aux déciles (1, 11, 21... 91)
            deciles = data_avg(1:10:end); 
            fprintf('  %-10s : %s\n', m_name, num2str(deciles, '%.2f '));
        end
    end

    % Enregistrement de la matrice de synergie globale pour cette condition
    if ~isempty(condition_matrix)
        SYNERGY_Adulte.(condition).matrix = condition_matrix;
        SYNERGY_Adulte.(condition).muscle_names = muscle_names;
    end
end

%% ===================== SAUVEGARDE FINALE ================================

suffix = upper(groupe_a_etudier);
fprintf('\nSauvegarde des fichiers consolidés (%s)...\n', suffix);

% Sauvegarde des 3 fichiers clés
save(['SPM-EMG-' suffix '.mat'], 'mean_cycles_combined');
save(['SYNERGY-' suffix '.mat'], 'SYNERGY_Adulte');
save(['PARTICIPANT-CYCLES-' suffix '.mat'], 'normalized_participant_cycles', 'participant_tracking');

fprintf('=== PROCESSUS DE FUSION TERMINÉ ===\n');