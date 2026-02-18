%% CALCUL EN BATCH DES MATRICES DE COHÉRENCE SIGNIFICATIVE (CYCLE ENTIER)
%
% OBJECTIF :
%   Ce script parcourt tous les fichiers Coherence_<PID>.mat d'un dossier
%   et calcule, pour chaque paire de muscles et chaque condition/côté, une
%   nouvelle matrice de cohérence ne conservant que la partie statistiquement
%   significative :
%
%       CoherenceSignif(f,t) = max( Coherence(f,t) - seuil, 0 )
%
%   Cette opération de rectification élimine les valeurs sous le seuil de
%   Rosenberg (calculé dans le script 7) sans modifier la matrice originale.
%   Le résultat est stocké dans un nouveau champ 'CoherenceSignif_<m1>_<m2>'
%   et le fichier .mat est mis à jour en place.
%
%   Seules les matrices de cycle complet sont traitées (champs de la forme
%   'Coherence_<m1>_<m2>' à exactement 3 parties). Les champs de sous-phases
%   ('Coherence_Phase_m1_m2') sont intentionnellement ignorés.
%
% ENTRÉES :
%   - Coherence_<PID>.mat : Fichiers de résultats issus du script 7,
%                           contenant la structure DATA avec, pour chaque
%                           condition/côté/paire :
%       * Coherence_<m1>_<m2>  : Matrice de cohérence [nFreq × 1000]
%       * Seuil_<m1>_<m2>      : Seuil de significativité scalaire
%
% SORTIES :
%   - Coherence_<PID>.mat : Même fichier, mis à jour avec les nouveaux champs :
%       * CoherenceSignif_<m1>_<m2> : Matrice de cohérence significative
%                                     [nFreq × 1000], valeurs ≥ 0
%       * meta.CoherenceSignif_FullCycle_Method : Description de la méthode
%       * meta.timestamp_CoherenceSignif_FullCycle : Horodatage du calcul
%
% OPTIONS :
%   - make_backup (bool) : Si true, crée une copie de sauvegarde du fichier
%                          avant modification (nommée *_BACKUP_<timestamp>)
%   - dry_run     (bool) : Si true, effectue tous les calculs sans écrire
%                          sur disque (mode test)
%
% NOTES :
%   - Les matrices originales (Coherence_*) ne sont jamais modifiées.
%   - Une erreur sur un fichier n'arrête pas le traitement des suivants.
%   - Un récapitulatif global (fichiers OK/ERREUR, matrices créées/sautées)
%     est affiché à la fin du batch.
% -------------------------------------------------------------------------

clear; clc; close all;

%% ===================== PARAMÈTRES =======================================

% Dossier contenant tous les fichiers Coherence_<PID>.mat à traiter
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';

% Options de traitement
make_backup = false;  % true : crée un fichier *_BACKUP avant modification
dry_run     = false;  % true : calcule mais ne sauvegarde pas (mode test)

assert(isfolder(all_dir), 'Dossier introuvable: %s', all_dir);

%% ===================== DÉTECTION DES FICHIERS ===========================

files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(files)
    fprintf('Aucun fichier Coherence_*.mat trouvé dans: %s\n', all_dir);
    return;
end
fprintf('>> %d fichier(s) détecté(s) dans %s\n', numel(files), all_dir);

%% ===================== INITIALISATIONS DES COMPTEURS GLOBAUX ============

global_created   = 0; % Nb de matrices CoherenceSignif créées sur l'ensemble du batch
global_skipped   = 0; % Nb de paires ignorées (seuil manquant ou matrice invalide)
global_files_ok  = 0; % Nb de fichiers traités avec succès
global_files_err = 0; % Nb de fichiers ayant généré une erreur
batch_notes      = {}; % Journal centralisé des avertissements

%% ===================== BOUCLE PRINCIPALE SUR LES FICHIERS ===============

for iFile = 1:numel(files)
    fpath = fullfile(files(iFile).folder, files(iFile).name);
    fprintf('\n================= [%d/%d] %s =================\n', ...
        iFile, numel(files), files(iFile).name);

    try
        % --- Chargement ---
        S = load(fpath);
        if ~isfield(S, 'DATA')
            error('Le fichier ne contient pas la variable DATA.');
        end
        DATA = S.DATA;

        % --- Filtrage des conditions valides ---
        condNames = fieldnames(DATA);
        validMask = false(size(condNames));
        for k = 1:numel(condNames)
            nm = condNames{k};
            validMask(k) = isstruct(DATA.(nm)) && ...
                           (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
        end
        condNames = condNames(validMask);
        sides     = {'left','right'};

        n_done = 0; n_skip = 0; local_notes = {};

        % --- Traitement par condition et côté ---
        for iC = 1:numel(condNames)
            condName = condNames{iC};

            % Annotation des métadonnées : méthode et horodatage du traitement
            try
                DATA.(condName).meta.CoherenceSignif_FullCycle_Method = ...
                    ['Element-wise NewC = max(Coherence - seuil, 0). ', ...
                     'Only full-cycle fields processed (Coherence_m1_m2).'];
                DATA.(condName).meta.timestamp_CoherenceSignif_FullCycle = ...
                    datestr(now,'yyyy-mm-dd HH:MM:SS');
            end

            for s = 1:numel(sides)
                sideStr = sides{s};
                if ~isfield(DATA.(condName), sideStr), continue; end
                Str = DATA.(condName).(sideStr);

                fns = fieldnames(Str);
                for iF = 1:numel(fns)
                    fn = fns{iF};
                    if ~startsWith(fn, 'Coherence'), continue; end

                    % --- Sélection du cycle complet uniquement ---
                    % Format attendu : 'Coherence_<m1>_<m2>' (exactement 3 parties)
                    % Les champs de sous-phases (4 parties) sont ignorés
                    parts = strsplit(fn, '_');
                    if numel(parts) ~= 3
                        continue; % Champ de sous-phase → ignoré intentionnellement
                    end
                    m1 = parts{2}; m2 = parts{3};

                    seuil_field = sprintf('Seuil_%s_%s', m1, m2);
                    outfield    = sprintf('CoherenceSignif_%s_%s', m1, m2);

                    % Vérification de la présence du seuil
                    if ~isfield(Str, seuil_field) || isempty(Str.(seuil_field))
                        n_skip = n_skip + 1;
                        local_notes{end+1} = sprintf('%s | %s | %s : seuil manquant (%s)', ...
                            condName, sideStr, fn, seuil_field);
                        continue;
                    end

                    % Vérification de la validité de la matrice de cohérence
                    if ~isfield(Str, fn) || isempty(Str.(fn)) || ndims(Str.(fn)) ~= 2
                        n_skip = n_skip + 1;
                        local_notes{end+1} = sprintf( ...
                            '%s | %s | %s : matrice absente/vide ou non 2D (size=%s)', ...
                            condName, sideStr, fn, mat2str(size(Str.(fn))));
                        continue;
                    end

                    C     = Str.(fn);
                    seuil = Str.(seuil_field);

                    % --- Calcul de la cohérence significative ---
                    % Soustraction élément par élément du seuil de Rosenberg,
                    % puis rectification : toute valeur négative → 0
                    % (cohérence non significative traitée comme nulle)
                    NewC = C - seuil;
                    NewC(NewC < 0) = 0;

                    % Écriture dans un NOUVEAU champ sans toucher à l'original
                    DATA.(condName).(sideStr).(outfield) = NewC;

                    n_done = n_done + 1;
                    fprintf('   ✓ %s → %s | seuil=%.4f | size=%s\n', ...
                        fn, outfield, seuil, mat2str(size(NewC)));
                end
            end
        end

        fprintf('   -> Nouvelles matrices créées: %d | Sautées: %d\n', n_done, n_skip);

        % --- Sauvegarde (avec option backup) ---
        if ~dry_run
            if make_backup
                [p,n,e] = fileparts(fpath);
                backup_file = fullfile(p, [n '_BACKUP_' datestr(now,'yyyymmdd_HHMMSS') e]);
                copyfile(fpath, backup_file);
                fprintf('   -> Backup créé: %s\n', backup_file);
            end
            save(fpath, 'DATA', '-v7.3'); % Mise à jour en place
            fprintf('   -> Sauvegardé: %s\n', fpath);
        else
            fprintf('   -> DRY RUN (aucune écriture disque)\n');
        end

        % Mise à jour des compteurs globaux
        global_created  = global_created  + n_done;
        global_skipped  = global_skipped  + n_skip;
        batch_notes     = [batch_notes, local_notes];
        global_files_ok = global_files_ok + 1;

    catch ME
        % Une erreur sur un fichier n'interrompt pas le batch
        global_files_err = global_files_err + 1;
        warnmsg = sprintf('ERREUR fichier %s: %s', files(iFile).name, ME.message);
        fprintf('   !! %s\n', warnmsg);
        batch_notes{end+1} = warnmsg;
    end
end

%% ===================== RÉCAPITULATIF GLOBAL =============================

fprintf('\n================= RÉCAP BATCH =================\n');
fprintf('Fichiers traités OK    : %d\n', global_files_ok);
fprintf('Fichiers en ERREUR     : %d\n', global_files_err);
fprintf('Matrices CoherenceSignif créées : %d\n', global_created);
fprintf('Paires ignorées        : %d\n', global_skipped);

if ~isempty(batch_notes)
    fprintf('\nNotes/Warnings (%d):\n', numel(batch_notes));
    % Affichage limité aux 30 premières notes pour ne pas noyer la console
    max_show = min(30, numel(batch_notes));
    for i = 1:max_show
        fprintf('  - %s\n', batch_notes{i});
    end
    if numel(batch_notes) > max_show
        fprintf('  ... (+%d supplémentaires)\n', numel(batch_notes) - max_show);
    end
end