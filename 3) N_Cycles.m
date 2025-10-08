%% AFFICHER LE NOMBRE DE CYCLES VALIDES PAR JAMBE, PAR CONDITION ET PAR PARTICIPANT
clc; clear; close all;

% Dossier contenant les fichiers .mat
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS');

% Récupérer tous les fichiers *_MATRIX.mat
files = dir('*_MATRIX.mat');

% --- Créer un tableau typé avec les bons noms de colonnes dès le départ
CycleTable = table('Size',[0 4], ...
    'VariableTypes', {'string','string','string','double'}, ...
    'VariableNames', {'Participant','Condition','Jambe','Nb_Cycles'});

% Parcourir chaque fichier
for i = 1:numel(files)
    filename = files(i).name;

    % Charger prudemment
    S = load(filename, 'CYCLES_COUNT');
    if ~isfield(S, 'CYCLES_COUNT')
        warning('Variable CYCLES_COUNT absente dans %s', filename);
        continue
    end
    CYCLES_COUNT = S.CYCLES_COUNT;

    % Extraire le nom du participant
    [~, name] = fileparts(filename);
    participant = erase(name, '_MATRIX');

    % Vérifier la structure
    if ~isfield(CYCLES_COUNT, participant)
        warning('Pas de champ %s dans CYCLES_COUNT de %s', participant, filename);
        continue
    end

    conditions = fieldnames(CYCLES_COUNT.(participant));
    for c = 1:numel(conditions)
        cond = conditions{c};
        jambes = fieldnames(CYCLES_COUNT.(participant).(cond));
        for j = 1:numel(jambes)
            leg = jambes{j};
            n_cycles = CYCLES_COUNT.(participant).(cond).(leg);

            % Si ce n’est pas un scalaire numérique (par ex. un vecteur d’indices), prendre sa taille
            if ~(isnumeric(n_cycles) && isscalar(n_cycles))
                n_cycles = numel(n_cycles);
            end

            % Ajouter une ligne (types respectés)
            CycleTable = [CycleTable; {string(participant), string(cond), string(leg), double(n_cycles)}];
        end
    end
end

% (Optionnel) trier pour lisibilité
CycleTable = sortrows(CycleTable, {'Participant','Condition','Jambe'});

% Afficher le tableau
disp(CycleTable);

% Chemin de sauvegarde
outdir = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results';
if ~exist(outdir, 'dir'); mkdir(outdir); end
output_path = fullfile(outdir, 'CycleTable.xlsx');

% Sauvegarde en Excel
writetable(CycleTable, output_path, 'FileType', 'spreadsheet');
fprintf('Fichier sauvegardé ici : %s\n', output_path);