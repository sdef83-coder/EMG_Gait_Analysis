%% EXTRACTION ET CONCATÉNATION DES CYCLES EMG BRUTS (.c3d -> .mat)
%
% OBJECTIF :
% Ce script automatise l'extraction des signaux EMG bruts à partir de fichiers 
% de capture de mouvement (.c3d). Il découpe le signal en cycles (Gait Cycles)
% basés sur les évènements de contact au sol (Heel Strikes) et concatène ces 
% cycles en un seul vecteur continu par muscle, condition et jambe.
%
% ENTRÉES :
%   - Fichiers .c3d : Contenant les données analogiques (EMG) et les événements.
%   - Association.m : Script de mapping liant les labels EMG aux capteurs réels.
%   - Fonctions : indiceLeft.m, indiceRight.m (détection des cycles).
%
% SORTIES :
%   - Fichiers .mat : Un fichier par participant contenant la structure 
%     'CYCLES_SIGNAL_RAW' (Signal concaténé et nombre de cycles).
%   - Images .png : Visualisation de contrôle du signal concaténé pour chaque muscle.
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CONFIGURATION & SÉLECTION ========================

% Mode de sélection des participants :
% 'auto'            : Détecte tous les .c3d dans DATA_DIR
% 'manual_name'     : Traite uniquement l'ID défini dans MANUAL_PID
% 'from_file_dialog': Ouvre une fenêtre de sélection de fichier
PARTICIPANT_MODE = 'manual_name';  
MANUAL_PID       = 'CTL_19'; 

% Chemins d'accès (À adapter selon l'environnement de travail)
DATA_DIR  = 'C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Data\.c3d\enfants';
OUT_DIR   = 'C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Data\.mat_raw_EMG_data';
PNG_DIR   = 'C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Data\.png_visualisation_EMG\Bruts';

% Paramètres expérimentaux
Essai       = {'01','02','03','04', '05','06','07','08','09','10'};
muscles     = {'EMG_TAprox','EMG_TAdist','EMG_SOL','EMG_GM','EMG_VL','EMG_RF','EMG_ST','EMG_GMED'};
Conditions  = {'Plat','Medium','High'};
jambes      = {'left','right'};
FreqVicon   = 100; % Fréquence d'acquisition vidéo (Hz)

%% ===================== PRÉPARATION DU WORKSPACE =========================

% Ajout des fonctions utilitaires au path Matlab
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Script\Functions'));

% Vérification et création des dossiers de sortie
if ~exist(DATA_DIR, 'dir'), error('DATA_DIR inexistant: %s', DATA_DIR); end
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
if ~exist(PNG_DIR, 'dir'), mkdir(PNG_DIR); end

% Identification des fichiers disponibles
allc3d = dir(fullfile(DATA_DIR, '*.c3d'));
if isempty(allc3d), error('Aucun fichier .c3d trouvé.'); end

% Pattern de nommage : PID_Condition_Essai.c3d
pid_re = '^(?<pid>[A-Za-z0-9]+_[A-Za-z0-9]+)_(?<cond>[A-Za-z]+)_(?<ess>\d{2})\.c3d$';

%% ===================== SÉLECTION DES PARTICIPANTS =======================

switch lower(PARTICIPANT_MODE)
    case 'auto'
        pid_set = {};
        for k = 1:numel(allc3d)
            tok = regexp(allc3d(k).name, pid_re, 'names');
            if ~isempty(tok), pid_set{end+1} = tok.pid; end
        end
        Participants = unique(pid_set);
    case 'manual_name'
        Participants = {MANUAL_PID};
    case 'from_file_dialog'
        [f, p] = uigetfile(fullfile(DATA_DIR, '*.c3d'), 'Sélectionner un fichier');
        tok = regexp(f, pid_re, 'names');
        Participants = {tok.pid};
end

fprintf('Participants à traiter (%d): %s\n', numel(Participants), strjoin(Participants, ', '));

%% ===================== BOUCLE DE TRAITEMENT PRINCIPALE ==================

for iP = 1:numel(Participants)
    pid = Participants{iP};
    fprintf('\n>>> Traitement : %s <<<\n', pid);

    % Chargement du mapping spécifique au participant (via Association.m)
    Participant = {pid}; 
    clear sensor_association_left sensor_association_right
    run Association.m               

    % Structure de données finale pour le participant
    CYCLES_SIGNAL_RAW = struct();
    CYCLES_SIGNAL_RAW.(pid) = struct();

    for iC = 1:numel(Conditions)
        cond = Conditions{iC};
        
        for j = 1:numel(jambes)
            leg = jambes{j};
            sensor_association = if(strcmp(leg,'left'), sensor_association_left, sensor_association_right);

            for m = 1:numel(muscles)
                mus = muscles{m};
                if ~isfield(sensor_association, mus), continue; end

                concat_sig = [];
                n_cycles_total = 0;

                % Parcours des essais pour extraire et concaténer les cycles
                for iE = 1:numel(Essai)
                    c3dname = sprintf('%s_%s_%s.c3d', pid, cond, Essai{iE});
                    fullpth = fullfile(DATA_DIR, c3dname);
                    if ~isfile(fullpth), continue; end

                    try
                        % Lecture des données C3D via BTK
                        acq     = btkReadAcquisition(fullpth);
                        analogs = btkGetAnalogs(acq);
                        FreqS   = btkGetAnalogFrequency(acq);
                        
                        sensor_name = sensor_association.(mus);
                        if ~isfield(analogs, sensor_name), continue; end

                        EMG_signal = analogs.(sensor_name);

                        % Identification des indices de début/fin de cycles (Heel Strikes)
                        if strcmp(leg,'left')
                            HS = indiceLeft(acq, analogs, FreqS, FreqVicon, EMG_signal);
                        else
                            HS = indiceRight(acq, analogs, FreqS, FreqVicon, EMG_signal);
                        end

                        % Découpe et accumulation du signal
                        for k = 1:(numel(HS)-1)
                            i_start = max(1, floor(HS(k)));
                            i_end   = min(numel(EMG_signal), ceil(HS(k+1)));
                            if i_end > i_start
                                concat_sig = [concat_sig; EMG_signal(i_start:i_end)];
                                n_cycles_total = n_cycles_total + 1;
                            end
                        end
                    catch ME
                        warning('Erreur sur %s : %s', c3dname, ME.message);
                    end
                end 

                % Enregistrement des données si des cycles ont été trouvés
                if ~isempty(concat_sig)
                    if ~isfield(CYCLES_SIGNAL_RAW.(pid), cond), CYCLES_SIGNAL_RAW.(pid).(cond) = struct(); end
                    if ~isfield(CYCLES_SIGNAL_RAW.(pid).(cond), leg), CYCLES_SIGNAL_RAW.(pid).(cond).(leg) = struct(); end

                    catField = [mus '_Concat'];
                    CYCLES_SIGNAL_RAW.(pid).(cond).(leg).(catField).signal   = concat_sig;
                    CYCLES_SIGNAL_RAW.(pid).(cond).(leg).(catField).n_cycles = n_cycles_total;

                    % Génération d'une figure de contrôle qualité
                    try
                        fig = figure('Visible','off','Color','w');
                        plot(concat_sig, 'LineWidth', 0.5); 
                        title(sprintf('EMG Brut Concaténé | %s | %s | %s', pid, leg, mus));
                        xlabel('Samples'); ylabel('Amplitude');
                        
                        saveas(fig, fullfile(PNG_DIR, sprintf('%s_%s_%s_%s.png', pid, cond, leg, mus)));
                        close(fig);
                    catch, end
                end
            end 
        end 
    end 

    % Sauvegarde physique du fichier .mat pour le participant
    save_file = fullfile(OUT_DIR, sprintf('%s_RAW_CYCLES_CONCAT.mat', pid));
    save(save_file, 'CYCLES_SIGNAL_RAW', '-v7.3');
    fprintf('Fichier sauvegardé: %s\n', save_file);
end

fprintf('\n=== FIN DE L''EXTRACTION ===\n');