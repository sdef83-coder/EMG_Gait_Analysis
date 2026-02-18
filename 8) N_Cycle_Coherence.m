%% EXTRACTION DU NOMBRE DE CYCLES UTILISÉS PAR ANALYSE DE COHÉRENCE
%
% OBJECTIF :
%   Ce script parcourt l'ensemble des fichiers de résultats de cohérence
%   (Coherence_<PID>.mat) générés par le script 7, et en extrait le nombre
%   de cycles effectivement utilisés pour chaque combinaison participant ×
%   condition × paire de muscles.
%
%   La valeur extraite est L_target_eff, soit le nombre de cycles après
%   égalisation inter-conditions (cf. script 7). Ce nombre est d'abord
%   recherché dans les métadonnées (DATA.(cond).meta), puis reconstruit
%   par sommation des côtés gauche et droit si les métadonnées sont absentes.
%
%   Le résultat est exporté dans un fichier CSV structuré, utilisable
%   directement pour les analyses statistiques ou la vérification de
%   l'équilibre du plan expérimental.
%
% ENTRÉES :
%   - Coherence_<PID>.mat : Fichiers de résultats issus du script 7,
%                           un par participant, contenant la structure DATA.
%                           Champ principal utilisé :
%                             DATA.(cond).meta.L_target_eff_<m1>_<m2>
%                           (ou DATA.(cond).left/right.L_<m1>_<m2> en fallback)
%
% SORTIES :
%   - N_Cycle_Coherence.csv : Tableau à 4 colonnes :
%       * Participant  : Identifiant du participant (ex : 'CTL_63')
%       * Condition    : Condition expérimentale ('Plat', 'Medium', 'High')
%       * Pair         : Paire de muscles (ex : 'TAprox_TAdist')
%       * L_target_eff : Nombre de cycles utilisés pour cette combinaison
%
% NOTES :
%   - Les fichiers .mat sont détectés automatiquement dans IN_DIR.
%   - Si L_target_eff est absent des métadonnées, le script tente une
%     reconstruction par sommation des champs L_<pair> sur les côtés
%     gauche et droit (comportement de secours / fallback).
%   - Les lignes avec valeur manquante (NaN) ne sont pas exportées.
%   - Le tableau final est trié par Participant, Condition, puis Pair.
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CHEMINS & VÉRIFICATIONS ==========================

% Dossier contenant tous les fichiers Coherence_<PID>.mat à traiter
IN_DIR  = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';

% Dossier de sortie pour le fichier CSV
OUT_DIR = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
OUT_CSV = fullfile(OUT_DIR, 'N_Cycle_Coherence.csv');

% Vérifications des dossiers
if ~exist(IN_DIR,'dir'),  error('Dossier introuvable: %s', IN_DIR); end
if ~exist(OUT_DIR,'dir'), mkdir(OUT_DIR); end

% Fonction anonyme pour extraire l'identifiant du participant depuis
% le nom de fichier (ex: 'Coherence_CTL_63.mat' → 'CTL_63')
extract_id = @(fn) string(regexprep(fn, '^Coherence_([^\.]+)\.mat$', '$1'));

% ===================== INITIALISATION DE LA TABLE DE SORTIE =============

% Table vide avec les types et noms de colonnes attendus
T = table('Size',[0 4], ...
          'VariableTypes', {'string','string','string','double'}, ...
          'VariableNames', {'Participant','Condition','Pair','L_target_eff'});

%% ===================== BOUCLE SUR LES FICHIERS .mat =====================

files = dir(fullfile(IN_DIR, 'Coherence_*.mat'));

for i = 1:numel(files)
    fpath = fullfile(files(i).folder, files(i).name);

    % Chargement du fichier — on ne lit que la variable 'DATA'
    S = load(fpath, 'DATA');
    if ~isfield(S,'DATA')
        warning('Variable DATA absente dans : %s', files(i).name);
        continue;
    end
    DATA = S.DATA;
    pid  = extract_id(files(i).name); % Identifiant du participant

    % --- Filtrage des conditions valides ---
    % Conserve uniquement les champs de DATA qui correspondent à des
    % conditions expérimentales (structs contenant 'left', 'right' ou 'meta')
    condNames = fieldnames(DATA);
    keep = false(size(condNames));
    for k = 1:numel(condNames)
        nm = condNames{k};
        keep(k) = isstruct(DATA.(nm)) && ...
                  (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right') || isfield(DATA.(nm),'meta'));
    end
    condNames = condNames(keep);

    for ic = 1:numel(condNames)
        cond = condNames{ic};

        % Récupération du sous-struct de métadonnées (peut être absent)
        meta_here = [];
        if isfield(DATA.(cond),'meta') && isstruct(DATA.(cond).meta)
            meta_here = DATA.(cond).meta;
        end

        % --- Détection des paires de muscles disponibles ---

        % Stratégie 1 (prioritaire) : liste les paires via les champs
        % 'L_target_eff_<m1>_<m2>' présents dans les métadonnées
        pair_list = {};
        if ~isempty(meta_here)
            mf   = fieldnames(meta_here);
            mask = startsWith(mf, 'L_target_eff_');
            keys = mf(mask);
            for j = 1:numel(keys)
                % Retire le préfixe pour ne garder que '<m1>_<m2>'
                pair_list{end+1} = extractAfter(keys{j}, 'L_target_eff_'); 
            end
            pair_list = unique(pair_list);
        end

        % Stratégie 2 (fallback) : déduit les paires depuis les champs
        % 'L_<m1>_<m2>' présents dans les sous-structs left/right
        if isempty(pair_list)
            tmp = {};
            for side = {'left','right'}
                s = side{1};
                if isfield(DATA.(cond), s)
                    fns = fieldnames(DATA.(cond).(s));
                    lm  = startsWith(fns, 'L_');
                    tmp = [tmp, extractAfter(fns(lm), 'L_')']; 
                end
            end
            pair_list = unique(tmp);
        end

        % Aucune paire détectée pour cette condition → passe à la suivante
        if isempty(pair_list), continue; end

        % --- Extraction de L_target_eff pour chaque paire ---
        for ip = 1:numel(pair_list)
            pair = pair_list{ip};
            val  = NaN;

            % Lecture directe depuis les métadonnées (cas nominal)
            fmeta = ['L_target_eff_' pair];
            if ~isempty(meta_here) && isfield(meta_here, fmeta)
                val = double(meta_here.(fmeta));

            else
                % Fallback : reconstruction par sommation gauche + droit
                % Si les deux côtés ont été traités séparément, on additionne
                % les L_<pair> de chaque côté pour obtenir le total
                fL = ['L_' pair];
                vL = NaN; vR = NaN;
                if isfield(DATA.(cond),'left')  && isfield(DATA.(cond).left,  fL)
                    vL = double(DATA.(cond).left.(fL));
                end
                if isfield(DATA.(cond),'right') && isfield(DATA.(cond).right, fL)
                    vR = double(DATA.(cond).right.(fL));
                end
                if ~isnan(vL) || ~isnan(vR)
                    val = nansum([vL, vR]); % Somme en ignorant les NaN individuels
                end
            end

            % Ajout à la table uniquement si une valeur valide a été trouvée
            if ~isnan(val)
                T = [T; {pid, string(cond), string(pair), val}]; %#ok<AGROW>
            end
        end
    end
end

%% ===================== EXPORT CSV =======================================

% Tri final par participant, condition, puis paire pour faciliter la lecture
T = sortrows(T, {'Participant','Condition','Pair'});
writetable(T, OUT_CSV);
fprintf('OK: %s\n', OUT_CSV);

%% ===================== FONCTIONS LOCALES ================================

function id = local_extract_id(fname)
% LOCAL_EXTRACT_ID  Extrait l'identifiant du participant depuis le nom de fichier.
%
%   Entrée  : fname  — nom du fichier (ex: 'Coherence_CTL_63.mat')
%   Sortie  : id     — identifiant string (ex: 'CTL_63'), ou "" si non trouvé
%
%   Utilisée en alternative à la fonction anonyme extract_id du script
%   principal (conservée ici pour compatibilité).

    id  = "";
    tok = regexp(fname, 'Coherence_(?<id>CTL_[^\.]+)\.mat$', 'names');
    if ~isempty(tok) && isfield(tok,'id')
        id = string(tok.id);
    end
end

function n = local_get_n_cycles(x)
% LOCAL_GET_N_CYCLES  Déduit le nombre de cycles depuis différents formats.
%
%   Cette fonction robuste accepte plusieurs types de données pouvant
%   représenter un nombre de cycles dans la structure DATA, et retourne
%   toujours un scalaire numérique (ou NaN si la déduction échoue).
%
%   Entrée : x — valeur à interpréter, de type :
%       - Scalaire numérique : retourné tel quel
%       - Vecteur numérique  : longueur du vecteur (numel)
%       - Logique scalaire   : converti en 0 ou 1
%       - Logique vecteur    : nombre de valeurs vraies (nnz)
%       - Cell               : nombre d'éléments (numel)
%       - Struct             : cherche un champ connu parmi une liste de
%                              candidats, puis tente le champ unique si struct
%                              à un seul champ
%
%   Sortie : n — nombre de cycles (double scalaire), ou NaN si indéterminé

    n = NaN;

    if isnumeric(x)
        if isscalar(x), n = double(x); else, n = numel(x); end
        return
    end

    if islogical(x)
        % Vecteur logique → nombre de cycles valides (vrais)
        if isscalar(x), n = double(x~=0); else, n = nnz(x); end
        return
    end

    if iscell(x)
        n = numel(x); return
    end

    if isstruct(x)
        % Liste des noms de champs usuels susceptibles de contenir N cycles
        candidates = {'Nb_Cycles','N','n','count','Valid_Cycles','valid_cycles', ...
                      'indices','idx','cycles','cycle_idx','valid_idx'};
        for c = 1:numel(candidates)
            f = candidates{c};
            if isfield(x, f)
                val = x.(f);
                if isnumeric(val) || islogical(val)
                    if isscalar(val), n = double(val); else, n = numel(val); end
                elseif iscell(val)
                    n = numel(val);
                end
                if ~isnan(n), return, end
            end
        end

        % Dernier recours : struct à un seul champ contenant une valeur numérique
        fns = fieldnames(x);
        if numel(fns) == 1
            val = x.(fns{1});
            if isnumeric(val) || islogical(val)
                if isscalar(val), n = double(val); else, n = numel(val); end
                return
            elseif iscell(val)
                n = numel(val); return
            end
        end
    end
end