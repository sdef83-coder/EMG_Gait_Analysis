%% === EXTRACTION CYCLES EMG BRUTS — CONCAT SEULEMENT =====================
% Sortie par participant (UN SEUL SIGNAL CONCATENE) :
%   CYCLES_SIGNAL_RAW.(pid).(Condition).(leg).([muscle '_Concat']).signal = [ ... ] (vector)

clc; clear; close all;

% ===================== CONFIG SELECTION (A ADAPTER) ======================
% Choisir comment définir les participants à traiter :
% 'auto'            : détecte tous les participants depuis les .c3d du dossier (attention à bien changer le nom du dossier DATA_DIR)
% 'manual_name'     : un seul participant, ID écrit ci-dessous (ex: 'CTL_07')
PARTICIPANT_MODE = 'manual_name';  % 'auto' | 'manual_name'
MANUAL_PID       = 'CTL_19';  % utilisé seulement si 'manual_name'

% ===================== PARAMETRES GENERAUX ===============================
DATA_DIR  = 'C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Data\.c3d\enfants';
OUT_DIR   = 'C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Data\.mat_raw_EMG_data';
PNG_DIR   = 'C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Data\.png_visualisation_EMG\Bruts';

Essai       = {'01','02','03','04', '05','06','07','08','09','10'}; %'05','06','07','08','09','10'
muscles     = {'EMG_TAprox','EMG_TAdist','EMG_SOL','EMG_GM','EMG_VL','EMG_RF','EMG_ST','EMG_GMED'};
Conditions  = {'Plat','Medium','High'};
jambes      = {'left','right'};
FreqVicon   = 100; 

% --- PATHS / DEPENDANCES -------------------------------------------------
addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Synergies_Projet_SurfaceIRR\Script\Functions'));   % Association.m, indiceLeft/Right, etc.
if ~exist(DATA_DIR, 'dir'), error('DATA_DIR inexistant: %s', DATA_DIR); end
if ~exist(OUT_DIR, 'dir'), mkdir(OUT_DIR); end
if ~exist(PNG_DIR, 'dir'), mkdir(PNG_DIR); end
cd(DATA_DIR);

% ===================== DETECTION DES FICHIERS ============================
allc3d = dir(fullfile(DATA_DIR, '*.c3d'));
if isempty(allc3d)
    error('Aucun fichier .c3d trouvé dans %s', DATA_DIR);
end

pid_re  = '^(?<pid>[A-Za-z0-9]+_[A-Za-z0-9]+)_(?<cond>[A-Za-z]+)_(?<ess>\d{2})\.c3d$';

% ===================== CHOIX DES PARTICIPANTS ============================
switch lower(PARTICIPANT_MODE)
    case 'auto'
        pid_set = {};
        for k = 1:numel(allc3d)
            tok = regexp(allc3d(k).name, pid_re, 'names');
            if ~isempty(tok)
                pid_set{end+1} = tok.pid; %#ok<AGROW>
            end
        end
        Participants = unique(pid_set);

        if isempty(Participants)
            error('Impossible d’extraire les IDs participants depuis les noms de fichiers. Vérifie le pattern.');
        end

    case 'manual_name'
        % Utilise l’ID fourni
        Participants = {MANUAL_PID};

        % Vérifie qu’au moins un .c3d existe pour cet ID
        has_file = any(contains({allc3d.name}, [MANUAL_PID '_']));
        if ~has_file
            warning('Aucun .c3d trouvé pour %s dans %s. Vérifie MANUAL_PID et DATA_DIR.', MANUAL_PID, DATA_DIR);
        end

    case 'from_file_dialog'
        % Laisse choisir un .c3d, et déduit l’ID
        [f, p] = uigetfile(fullfile(DATA_DIR, '*.c3d'), 'Sélectionne un .c3d pour inférer le participant');
        if isequal(f,0)
            error('Aucun fichier sélectionné.');
        end
        tok = regexp(f, pid_re, 'names');
        if isempty(tok)
            error('Le nom du fichier sélectionné ne matche pas le pattern attendu: %s', f);
        end
        Participants = {tok.pid};

    otherwise
        error('PARTICIPANT_MODE invalide: %s', PARTICIPANT_MODE);
end

fprintf('Participants à traiter (%d): %s\n', numel(Participants), strjoin(Participants, ', '));
fprintf('\n=== DÉBUT EXTRACTION ===\n');

% ===================== BOUCLE PRINCIPALE ================================
for iP = 1:numel(Participants)
    pid = Participants{iP};
    fprintf('\n--- Participant: %s ---\n', pid);

    % IMPORTANT: préparer Association.m (il attend Participant{1})
    Participant = {pid};            %#ok<NASGU>
    clear sensor_association_left sensor_association_right
    run Association.m               % doit définir sensor_association_left/right

    if ~exist('sensor_association_left','var') || ~exist('sensor_association_right','var')
        error('Association.m n’a pas créé sensor_association_left/right pour %s.', pid);
    end

    % Structure sortie pour ce participant
    CYCLES_SIGNAL_RAW = struct();
    CYCLES_SIGNAL_RAW.(pid) = struct();

    for iC = 1:numel(Conditions)
        cond = Conditions{iC};
        fprintf('  Condition: %s\n', cond);

        for j = 1:numel(jambes)
            leg = jambes{j};

            % Sélection de l’association selon la jambe
            if strcmp(leg,'left')
                sensor_association = sensor_association_left;
            else
                sensor_association = sensor_association_right;
            end

            for m = 1:numel(muscles)
                mus = muscles{m};

                if ~isfield(sensor_association, mus)
                    warning('    [SKIP] Pas d’association pour %s (%s, %s, %s).', mus, pid, cond, leg);
                    continue;
                end

                % === Un seul signal concaténé (tous cycles) ===
                concat_sig = [];
                n_cycles_total = 0;

                for iE = 1:numel(Essai)
                    ess = Essai{iE};
                    c3dname = sprintf('%s_%s_%s.c3d', pid, cond, ess);
                    fullpth = fullfile(DATA_DIR, c3dname);
                    if ~isfile(fullpth)
                        continue;
                    end

                    try
                        acq     = btkReadAcquisition(fullpth);
                        analogs = btkGetAnalogs(acq);
                        FreqS   = btkGetAnalogFrequency(acq);

                        sensor_name = sensor_association.(mus);
                        if ~isfield(analogs, sensor_name)
                            warning('    [SKIP] %s absent dans %s', sensor_name, c3dname);
                            continue;
                        end

                        EMG_signal = analogs.(sensor_name);  % *** BRUT ***

                        % Détection HS pour découpe des cycles
                        if strcmp(leg,'left')
                            HS = indiceLeft(acq, analogs, FreqS, FreqVicon, EMG_signal);
                        else
                            HS = indiceRight(acq, analogs, FreqS, FreqVicon, EMG_signal);
                        end

                        if numel(HS) < 2
                            continue;
                        end

                        % Concaténer chaque cycle bout-à-bout (sans NaN)
                        for k = 1:(numel(HS)-1)
                            i_start = max(1, floor(HS(k)));
                            i_end   = min(numel(EMG_signal), ceil(HS(k+1)));
                            if i_end > i_start
                                concat_sig = [concat_sig; EMG_signal(i_start:i_end)]; %#ok<AGROW>
                                n_cycles_total = n_cycles_total + 1;
                            end
                        end

                    catch ME
                        warning('    [ERR] %s -> %s', c3dname, ME.message);
                        continue;
                    end
                end % essais

                % Écriture si des cycles existent
                if ~isempty(concat_sig)
                    if ~isfield(CYCLES_SIGNAL_RAW.(pid), cond), CYCLES_SIGNAL_RAW.(pid).(cond) = struct(); end
                    if ~isfield(CYCLES_SIGNAL_RAW.(pid).(cond), leg), CYCLES_SIGNAL_RAW.(pid).(cond).(leg) = struct(); end

                    % Stocke UNIQUEMENT le signal concaténé
                    catField = [mus '_Concat'];
                    CYCLES_SIGNAL_RAW.(pid).(cond).(leg).(catField) = struct();
                    CYCLES_SIGNAL_RAW.(pid).(cond).(leg).(catField).signal      = concat_sig;
                    CYCLES_SIGNAL_RAW.(pid).(cond).(leg).(catField).n_cycles    = n_cycles_total;

                    fprintf('    %-10s (%s): %3d cycles concaténés (%d samples)\n', mus, leg, n_cycles_total, numel(concat_sig));

                    % === Sauvegarde de la figure (invisible) ===
                    try
                        fig = figure('Visible','off','Color','w','Name',sprintf('%s|%s|%s|%s', pid, cond, leg, mus));
                        ax  = axes(fig); %#ok<LAXES>
                        plot(ax, concat_sig, 'LineWidth', 0.5); grid(ax, 'on');
                        title(ax, sprintf('EMG brut concaténé — %s | %s | %s | %s', pid, cond, leg, strrep(mus,'_','\_')));
                        xlabel(ax, 'Samples'); ylabel(ax, 'Amplitude (brut)');

                        % nom de fichier PNG
                        png_name = sprintf('%s_%s_%s_%s_concat.png', pid, cond, leg, mus);
                        png_path = fullfile(PNG_DIR, png_name);

                        exportgraphics(ax, png_path, 'Resolution', 300);
                        close(fig);
                    catch MEpng
                        warning('    [PNG] Echec sauvegarde figure (%s | %s | %s | %s): %s', pid, cond, leg, mus, MEpng.message);
                    end
                else
                    fprintf('    %-10s (%s):   0 cycle\n', mus, leg);
                end

            end % muscles
        end % jambes
    end % conditions

    % Sauvegarde : 1 fichier par participant (concat uniquement)
    save_file = fullfile(OUT_DIR, sprintf('%s_RAW_CYCLES_CONCAT.mat', pid));
    save(save_file, 'CYCLES_SIGNAL_RAW', '-v7.3');
    fprintf('>>> Sauvegardé: %s\n', save_file);
end

fprintf('\n=== FIN EXTRACTION ===\n');