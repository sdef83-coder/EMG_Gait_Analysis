%% === BATCH: COHERENCE "SIGNIFICATIVE" (Cycle entier) POUR TOUS LES .MAT DU DOSSIER ===
% Parcourt tous les fichiers Coherence_*.mat dans ...\Coherence\ALL
% Pour chaque fichier:
%   - Charge DATA
%   - Pour chaque condition/côté:
%       * repère les champs Coherence_<m1>_<m2> (cycle complet uniquement)
%       * lit Seuil_<m1>_<m2> (scalaire)
%       * calcule NewC = max(Coherence - seuil, 0) élément par élément
%       * écrit dans CoherenceSignif_<m1>_<m2> (même taille que l’originale)
%   - Sauvegarde le .mat (après avoir créé un backup)
%
% Remarques:
%   - Les matrices d’origine ne sont pas modifiées.
%   - Les sous-phases (LoadingResponse/MidStance/PreSwing/Swing) sont ignorées.
%   - Un fichier en erreur n arrête pas le batch.

clear; clc; close all;

% ===== Paramètres dossier =====
all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';

% Options
make_backup = false;   % mettre true si j'en veux un
dry_run     = false;  % si true: ne sauvegarde pas, juste un "test"

assert(isfolder(all_dir), 'Dossier introuvable: %s', all_dir);

files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(files)
    fprintf('Aucun fichier Coherence_*.mat trouvé dans: %s\n', all_dir);
    return;
end

fprintf('>> %d fichier(s) détecté(s) dans %s\n', numel(files), all_dir);

% ===== Stats globales =====
global_created = 0;
global_skipped = 0;
global_files_ok = 0;
global_files_err = 0;

batch_notes = {};

% ===== Boucle fichiers =====
for iFile = 1:numel(files)
    fpath = fullfile(files(iFile).folder, files(iFile).name);
    fprintf('\n================= [%d/%d] %s =================\n', iFile, numel(files), files(iFile).name);

    try
        % Charger
        S = load(fpath);
        if ~isfield(S, 'DATA')
            error('Le fichier ne contient pas la variable DATA.');
        end
        DATA = S.DATA;

        % Déduire conditions
        condNames = fieldnames(DATA);
        validMask = false(size(condNames));
        for k = 1:numel(condNames)
            nm = condNames{k};
            validMask(k) = isstruct(DATA.(nm)) && (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
        end
        condNames = condNames(validMask);
        sides = {'left','right'};

        n_done = 0; n_skip = 0; local_notes = {};

        % Traiter chaque condition/côté
        for iC = 1:numel(condNames)
            condName = condNames{iC};

            % Meta
            try
                DATA.(condName).meta.CoherenceSignif_FullCycle_Method = ...
                    ['Element-wise NewC = max(Coherence - seuil, 0). ', ...
                     'Only full-cycle fields processed (Coherence_m1_m2).'];
                DATA.(condName).meta.timestamp_CoherenceSignif_FullCycle = datestr(now,'yyyy-mm-dd HH:MM:SS');
            end

            for s = 1:numel(sides)
                sideStr = sides{s};
                if ~isfield(DATA.(condName), sideStr), continue; end
                Str = DATA.(condName).(sideStr);

                fns = fieldnames(Str);
                for iF = 1:numel(fns)
                    fn = fns{iF};
                    if ~startsWith(fn, 'Coherence'), continue; end

                    % Cycle entier uniquement: "Coherence_<m1>_<m2>"
                    parts = strsplit(fn, '_');
                    if numel(parts) ~= 3
                        % Probable sous-phase -> on ignore
                        continue;
                    end
                    m1 = parts{2}; m2 = parts{3};

                    seuil_field = sprintf('Seuil_%s_%s', m1, m2);
                    outfield    = sprintf('CoherenceSignif_%s_%s', m1, m2);

                    if ~isfield(Str, seuil_field) || isempty(Str.(seuil_field))
                        n_skip = n_skip + 1;
                        local_notes{end+1} = sprintf('%s | %s | %s : seuil manquant (%s)', ...
                            condName, sideStr, fn, seuil_field);
                        continue;
                    end
                    if ~isfield(Str, fn) || isempty(Str.(fn)) || ndims(Str.(fn)) ~= 2
                        n_skip = n_skip + 1;
                        local_notes{end+1} = sprintf('%s | %s | %s : matrice absente/vide ou non 2D (size=%s)', ...
                            condName, sideStr, fn, mat2str(size(Str.(fn))));
                        continue;
                    end

                    C = Str.(fn);
                    seuil = Str.(seuil_field);

                    % Calcul élément par élément
                    NewC = C - seuil;
                    NewC(NewC < 0) = 0;

                    % Écrire nouvelle matrice (sans toucher à l'originale)
                    DATA.(condName).(sideStr).(outfield) = NewC;

                    n_done = n_done + 1;
                    fprintf('   ✓ %s → %s | seuil=%.4f | size=%s\n', fn, outfield, seuil, mat2str(size(NewC)));
                end
            end
        end

        fprintf('   -> Nouvelles matrices (créées ici): %d | Sautées: %d\n', n_done, n_skip);

        % Sauvegarde + backup
        if ~dry_run
            if make_backup
                [p,n,e] = fileparts(fpath);
                backup_file = fullfile(p, [n '_BACKUP_' datestr(now,'yyyymmdd_HHMMSS') e]);
                copyfile(fpath, backup_file);
                fprintf('   -> Backup créé: %s\n', backup_file);
            end
            save(fpath, 'DATA', '-v7.3');
            fprintf('   -> Sauvegardé: %s\n', fpath);
        else
            fprintf('   -> DRY RUN (aucune écriture)\n');
        end

        % Stats globales
        global_created = global_created + n_done;
        global_skipped = global_skipped + n_skip;
        batch_notes = [batch_notes, local_notes];
        global_files_ok = global_files_ok + 1;

    catch ME
        global_files_err = global_files_err + 1;
        warnmsg = sprintf('ERREUR fichier %s: %s', files(iFile).name, ME.message);
        fprintf('   !! %s\n', warnmsg);
        batch_notes{end+1} = warnmsg;
        % continue vers fichier suivant
    end
end

% ===== Récapitulatif =====
fprintf('\n================= RÉCAP BATCH =================\n');
fprintf('Fichiers OK:    %d\n', global_files_ok);
fprintf('Fichiers ERREUR:%d\n', global_files_err);
fprintf('Matrices créées:%d\n', global_created);
fprintf('Entrées sautées:%d\n', global_skipped);

if ~isempty(batch_notes)
    fprintf('\nNotes/Warnings (%d):\n', numel(batch_notes));
    % Affiche seulement les 30 premières pour ne pas noyer la console
    max_show = min(30, numel(batch_notes));
    for i = 1:max_show
        fprintf('  - %s\n', batch_notes{i});
    end
    if numel(batch_notes) > max_show
        fprintf('  ... (+%d supplémentaires)\n', numel(batch_notes) - max_show);
    end
end