% MAPPING EMG SUR TOUT LES MUSCLES SUR LE 1ER ET 10EME ESSAI DE CHAQUE SURFACE
% Penser à changer le groupe d'intérêt à la L6

clc; clear; close all;

% Dossier de base
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Data\\jeunes_enfants\');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Functions'));
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire'));

% Paramètres
Participant = {'CTL_63'};
Condition = {'Plat', 'Medium', 'High'};
Essais = {'01', '04'};
muscles = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'}; %'EMG_GMED'
jambes = {'left', 'right'};
 
run Association.m

FreqS = 2000;  % Hz

% Dossier de sauvegarde des figures
baseSavePath = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Fig\Activation';

for iP = 1:length(Participant)
    participant = Participant{iP};

    % Créer un dossier spécifique au participant s’il n’existe pas
    participantFolder = fullfile(baseSavePath, participant);
    if ~exist(participantFolder, 'dir')
        mkdir(participantFolder);
    end

    for iC = 1:length(Condition)
        for iM = 1:length(muscles)
            for iJ = 1:length(jambes)

                muscle = muscles{iM};
                jambe = jambes{iJ};

                if strcmp(jambe, 'left')
                    capteur = sensor_association_left.(muscle);
                else
                    capteur = sensor_association_right.(muscle);
                end

                signal_filtre_essai = cell(1,2);

                for iE = 1:2
                    essai = Essais{iE};
                    nom_fichier = sprintf('%s_%s_%s.c3d', participant, Condition{iC}, essai);

                    if isfile(nom_fichier)
                        data = btkReadAcquisition(nom_fichier);
                        analogs = btkGetAnalogs(data);
                        signal_brut = analogs.(capteur);
                        signal_filtre_essai{iE} = filtrage(signal_brut, FreqS, 20, 450);
                    else
                        warning('Fichier non trouvé : %s', nom_fichier);
                        signal_filtre_essai{iE} = [];
                    end
                end

                % Affichage figure
                fig = figure('Visible', 'off', ...
                    'Name', sprintf('%s - %s - %s - %s', participant, Condition{iC}, muscle, jambe), ...
                    'NumberTitle', 'off');

                for iE = 1:2
                    subplot(2,1,iE);
                    if ~isempty(signal_filtre_essai{iE})
                        plot(signal_filtre_essai{iE}, 'LineWidth', 1);
                        title(sprintf('Essai %s - %s - %s - %s', Essais{iE}, Condition{iC}, muscle, jambe), 'Interpreter', 'none');
                        xlabel('Temps (échantillons)');
                        ylabel('Amplitude EMG (a.u.)');
                        grid on;
                    else
                        title(sprintf('Essai %s manquant', Essais{iE}));
                    end
                end

                % Enregistrement figure
                figName = sprintf('%s_%s_%s_%s.png', participant, Condition{iC}, muscle, jambe);
                saveasPath = fullfile(participantFolder, figName);
                exportgraphics(fig, saveasPath, 'Resolution', 300);

                close(fig);  % Fermer la figure après sauvegarde pour libérer la mémoire
            end
        end
    end
end