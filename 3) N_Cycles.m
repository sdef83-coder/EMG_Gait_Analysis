%% RÉCAPITULATIF DU NOMBRE DE CYCLES VALIDES (EXPORT EXCEL)
%
% OBJECTIF :
% Ce script parcourt les résultats issus de la normalisation pour extraire le 
% nombre de cycles valides (après exclusion des outliers) pour chaque 
% participant, condition et jambe. Il génère un tableau récapitulatif 
% affiché dans la console et sauvegardé au format Excel (.xlsx).
%
% ENTRÉES :
%   - Fichiers *_MATRIX.mat : Contenant la structure 'CYCLES_COUNT'.
%
% SORTIES :
%   - CycleTable : Tableau Matlab (Table) affiché dans la console.
%   - CycleTable.xlsx : Fichier Excel exporté dans le dossier Results.
%
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CONFIGURATION DES CHEMINS ========================

% Dossier contenant les matrices de résultats
DATA_PATH = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results\Matrix\ORIGINALS';
% Dossier de destination pour l'export Excel
OUT_PATH  = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces_Irregulieres\Datas\Script\ActivationMusculaire\Results';

if ~exist(DATA_PATH, 'dir'), error('Dossier de données introuvable : %s', DATA_PATH); end
cd(DATA_PATH);

% Récupération de la liste des fichiers de résultats
files = dir('*_MATRIX.mat');
if isempty(files), error('Aucun fichier *_MATRIX.mat trouvé.'); end

% ===================== INITIALISATION DU TABLEAU ========================

% Création d'un tableau vide avec des types de colonnes définis pour la performance
% Colonnes : Nom Participant, Condition Expérimentale, Jambe, Nombre de Cycles
CycleTable = table('Size', [0 4], ...
    'VariableTypes', {'string', 'string', 'string', 'double'}, ...
    'VariableNames', {'Participant', 'Condition', 'Jambe', 'Nb_Cycles'});

%% ===================== EXTRACTION DES DONNÉES ===========================

fprintf('Analyse des fichiers en cours...\n');

for i = 1:numel(files)
    filename = files(i).name;

    % Chargement sélectif de la variable 'CYCLES_COUNT' uniquement
    S = load(filename, 'CYCLES_COUNT');
    if ~isfield(S, 'CYCLES_COUNT')
        warning('Variable CYCLES_COUNT absente dans %s. Fichier ignoré.', filename);
        continue
    end
    CYCLES_COUNT = S.CYCLES_COUNT;

    % Extraction de l'ID participant (ex: 'CTL_63' à partir de 'CTL_63_MATRIX.mat')
    [~, name] = fileparts(filename);
    participant = erase(name, '_MATRIX');

    % Vérification de la cohérence de la structure interne
    if ~isfield(CYCLES_COUNT, participant)
        warning('L''ID participant "%s" ne correspond pas au contenu du fichier %s', participant, filename);
        continue
    end

    % --- Boucles imbriquées pour parcourir la hiérarchie des données ---
    conditions = fieldnames(CYCLES_COUNT.(participant));
    for c = 1:numel(conditions)
        cond = conditions{c};
        
        jambes = fieldnames(CYCLES_COUNT.(participant).(cond));
        for j = 1:numel(jambes)
            leg = jambes{j};
            
            % Récupération de la valeur du nombre de cycles
            n_cycles = CYCLES_COUNT.(participant).(cond).(leg);

            % Sécurité : Si n_cycles est une liste d'indices au lieu d'un nombre, 
            % on calcule sa longueur pour obtenir le compte.
            if ~(isnumeric(n_cycles) && isscalar(n_cycles))
                n_cycles = numel(n_cycles);
            end

            % Ajout des données au tableau final
            newRow = {string(participant), string(cond), string(leg), double(n_cycles)};
            CycleTable = [CycleTable; newRow]; 
        end
    end
end

% ===================== TRI ET AFFICHAGE =================================

% Tri pour une lecture plus aisée (Participant > Condition > Jambe)
CycleTable = sortrows(CycleTable, {'Participant', 'Condition', 'Jambe'});

% Affichage du résultat final dans la fenêtre de commande
fprintf('\n--- RÉCAPITULATIF DES CYCLES VALIDES ---\n');
disp(CycleTable);

% ===================== EXPORTATION DES DONNÉES ==========================

if ~exist(OUT_PATH, 'dir'), mkdir(OUT_PATH); end
output_file = fullfile(OUT_PATH, 'CycleTable.xlsx');

% Sauvegarde au format Excel
try
    writetable(CycleTable, output_file, 'FileType', 'spreadsheet');
    fprintf('Fichier Excel généré avec succès :\n%s\n', output_file);
catch ME
    warning('Erreur lors de l''écriture du fichier Excel : %s', ME.message);
end

fprintf('\n=== OPÉRATION TERMINÉE ===\n');