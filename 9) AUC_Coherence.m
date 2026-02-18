%% CALCUL DE L'AIRE SIGNIFICATIVE SOUS LA COURBE DE COHÉRENCE (AUC)
%
% OBJECTIF :
%   Ce script calcule l'aire sous la courbe (AUC) de cohérence EMG-EMG pour
%   chaque combinaison condition × jambe × paire de muscles × sous-phase ×
%   bande fréquentielle. Seule la portion statistiquement significative de
%   la cohérence est intégrée : le seuil de significativité (calculé dans
%   le script 7) est soustrait au spectre de cohérence moyen, et les valeurs
%   négatives (non significatives) sont ramenées à zéro avant intégration.
%
%   L'intégration est réalisée par la méthode des trapèzes (trapz) sur
%   chaque bande fréquentielle d'intérêt (Alpha, Beta, Gamma).
%
%   Les résultats sont ajoutés directement dans la structure DATA du fichier
%   .mat existant, puis celui-ci est ré-enregistré (écrasement en place).
%
% ENTRÉES :
%   - Coherence_<PID>.mat : Fichier de résultats issu du script 7, contenant
%                           la structure DATA avec, pour chaque condition,
%                           jambe et paire de muscles :
%       * Coherence_<m1>_<m2>             : Cohérence sur le cycle complet
%       * Coherence_<Phase>_<m1>_<m2>     : Cohérence par sous-phase
%       * Seuil_<m1>_<m2>                 : Seuil de significativité (Rosenberg)
%       * Freq                            : Axe fréquentiel (Hz)
%
% SORTIES :
%   - Coherence_<PID>.mat : Même fichier, mis à jour avec les nouveaux champs :
%       * CoherenceArea_<Bande>_<m1>_<m2>              : AUC cycle complet
%       * CoherenceArea_<Phase>_<Bande>_<m1>_<m2>      : AUC par sous-phase
%         où Bande ∈ {Alpha, Beta, Gamma}
%              Phase ∈ {LoadingResponse, MidStance, PreSwing, Swing}
%       * meta.CoherenceArea_FreqBands   : Définition des bandes utilisées
%       * meta.CoherenceArea_Method      : Description de la méthode de calcul
%
% MÉTHODE DE CALCUL DE L'AUC SIGNIFICATIVE :
%   1. Moyenne temporelle de la cohérence sur la sous-phase : CohMean(f)
%   2. Soustraction du seuil : CohSignificant(f) = CohMean(f) - seuil
%   3. Mise à zéro des valeurs négatives (cohérence sous le seuil)
%   4. Remplissage des NaN éventuels (interpolation linéaire puis nearest)
%   5. Intégration numérique : AUC = trapz(f, CohSignificant(f)) sur la bande
%
% NOTES :
%   - Le fichier est lu ET réécrit au même emplacement (mise à jour en place).
%   - Les paires dont le seuil est absent sont ignorées et signalées.
%   - Les fréquences sont triées par ordre croissant avant intégration.
% -------------------------------------------------------------------------

clear; close all; clc;

%% ===================== PARAMÈTRES UTILISATEUR ===========================

% Identifiant du participant à traiter
participant_id = 'CTL_63'; % <- adapter si besoin

% Chemin racine des résultats de cohérence
coh_root = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(coh_root, 'ALL');

% Chemin du fichier à lire et à mettre à jour (même fichier en entrée et sortie)
infile  = fullfile(all_dir, ['Coherence_' participant_id '.mat']);
outfile = infile; % Réécriture en place

% Définition des bandes fréquentielles d'intérêt (Hz)
FreqBands = struct('Alpha',[8 12], 'Beta',[13 30], 'Gamma',[31 60]);

%% ===================== CHARGEMENT DU FICHIER ============================

if ~isfile(infile)
    error('Fichier introuvable: %s', infile);
end

L = load(infile);
if ~isfield(L, 'DATA')
    error('Le fichier ne contient pas la variable DATA: %s', infile);
end
DATA = L.DATA;

fprintf('>> Chargé: %s\n', infile);

%% ===================== FILTRAGE DES CONDITIONS VALIDES ==================

% Conserve uniquement les champs de DATA qui correspondent à des conditions
% expérimentales (structs contenant au moins un côté 'left' ou 'right')
condNames = fieldnames(DATA);
keep = false(size(condNames));
for k = 1:numel(condNames)
    nm = condNames{k};
    keep(k) = isstruct(DATA.(nm)) && ...
              (isfield(DATA.(nm),'left') || isfield(DATA.(nm),'right'));
end
condNames = condNames(keep);

%% ===================== BOUCLES PRINCIPALES ==============================

sides     = {'left','right'};
n_written = 0; % Compteur de champs AUC écrits avec succès
n_skipped = 0; % Compteur de paires ignorées (seuil manquant)
notes     = {}; % Journal des avertissements non bloquants

for iC = 1:numel(condNames)
    condName = condNames{iC};
    fprintf('\n== Condition: %s ==\n', condName);

    % Annotation des métadonnées de la condition avec la méthode utilisée
    % (utile pour la traçabilité lors des analyses ultérieures)
    try
        DATA.(condName).meta.CoherenceArea_FreqBands = FreqBands;
        DATA.(condName).meta.CoherenceArea_Method = ...
            'Mean over time (within phase), subtract significance threshold, set negative to 0, then trapz(Freq) over band. Freq sorted asc; NaNs filled.';
    end

    for s = 1:numel(sides)
        sideStr = sides{s};
        if ~isfield(DATA.(condName), sideStr), continue; end
        fprintf('-- Side: %s --\n', sideStr);
        Sside = DATA.(condName).(sideStr); % Raccourci vers le sous-struct du côté

        % Vérification de la présence de l'axe fréquentiel (indispensable)
        if ~isfield(Sside, 'Freq') || isempty(Sside.Freq)
            notes{end+1} = sprintf('%s | %s : Freq manquant', condName, sideStr);
            continue;
        end
        Freq = Sside.Freq(:); % Assure un vecteur colonne

        % --- Parcours de tous les champs de cohérence du côté ---
        % Deux formats de noms sont attendus :
        %   - 'Coherence_<m1>_<m2>'           → cycle complet (phase = '')
        %   - 'Coherence_<Phase>_<m1>_<m2>'   → sous-phase spécifique
        fns = fieldnames(Sside);
        for iF = 1:numel(fns)
            fn = fns{iF};
            if ~startsWith(fn,'Coherence'), continue; end

            % --- Décomposition du nom de champ ---
            parts = strsplit(fn,'_');
            phase = ''; m1 = ''; m2 = '';

            if numel(parts) == 3
                % Format : Coherence_<m1>_<m2> (cycle complet)
                m1 = parts{2}; m2 = parts{3};

            elseif numel(parts) == 4 && ...
                   any(strcmp(parts{2}, {'LoadingResponse','MidStance','PreSwing','Swing'}))
                % Format : Coherence_<Phase>_<m1>_<m2>
                phase = parts{2}; m1 = parts{3}; m2 = parts{4};
            else
                continue; % Champ de cohérence non reconnu → ignoré
            end

            % --- Récupération du seuil de significativité ---
            % Le seuil est stocké sous la forme 'Seuil_<m1>_<m2>'
            % Il est identique pour toutes les sous-phases d'une même paire
            seuil_field = sprintf('Seuil_%s_%s', m1, m2);

            if ~isfield(Sside, seuil_field) || isempty(Sside.(seuil_field))
                notes{end+1} = sprintf('%s | %s | %s : Seuil manquant (%s)', ...
                    condName, sideStr, fn, seuil_field);
                n_skipped = n_skipped + 1;
                continue;
            end
            seuil = Sside.(seuil_field);

            % Vérification de la cohérence dimensionnelle
            C = Sside.(fn); % Matrice [nFreq × nTemps]
            if isempty(C) || size(C,1) ~= numel(Freq)
                notes{end+1} = sprintf('%s | %s | %s : dims mismatch (C=%s, nFreq=%d)', ...
                    condName, sideStr, fn, mat2str(size(C)), numel(Freq));
                continue;
            end

            % --- Calcul de l'AUC significative ---

            % Étape 1 : Moyenne temporelle → spectre de cohérence moyen C(f)
            CohMean = mean(C, 2, 'omitnan');

            % Étape 2 : Soustraction du seuil (seuil de Rosenberg)
            % → isole la cohérence au-dessus du niveau de chance
            CohSignificant = CohMean - seuil;

            % Étape 3 : Rectification — les valeurs négatives (non significatives)
            % sont ramenées à zéro pour ne pas réduire l'AUC
            CohSignificant(CohSignificant < 0) = 0;

            fprintf('   %s: Seuil=%.4f, Orig max=%.4f, Signif max=%.4f\n', ...
                fn, seuil, max(CohMean), max(CohSignificant));

            % Étape 4 : Préparation de l'intégration
            % Tri par fréquence croissante (requis par trapz)
            [Fsort, idx] = sort(Freq, 'ascend');
            Csort = CohSignificant(idx);

            % Remplissage des NaN éventuels (interpolation linéaire, puis
            % nearest pour les extrémités)
            Csort = fillmissing(Csort, 'linear');
            Csort = fillmissing(Csort, 'nearest');

            % --- Étape 5 : Intégration par bande fréquentielle ---
            bandNames = fieldnames(FreqBands);
            for ib = 1:numel(bandNames)
                bname  = bandNames{ib};
                brange = FreqBands.(bname);

                % Sélection des fréquences appartenant à la bande
                mask     = (Fsort >= brange(1)) & (Fsort <= brange(2));
                area_val = NaN;

                if any(mask)
                    % Intégration numérique par la méthode des trapèzes
                    area_val = trapz(Fsort(mask), Csort(mask));
                end

                % Construction du nom du champ de sortie selon le format :
                %   Cycle complet : CoherenceArea_<Bande>_<m1>_<m2>
                %   Sous-phase    : CoherenceArea_<Phase>_<Bande>_<m1>_<m2>
                if isempty(phase)
                    outfield = sprintf('CoherenceArea_%s_%s_%s', bname, m1, m2);
                else
                    outfield = sprintf('CoherenceArea_%s_%s_%s_%s', phase, bname, m1, m2);
                end

                % Écriture de l'AUC dans la structure DATA
                DATA.(condName).(sideStr).(outfield) = area_val;
                n_written = n_written + 1;
            end
        end % Champs de cohérence
    end % Côtés
end % Conditions

%% ===================== RAPPORT DE TRAITEMENT ============================

fprintf('\n>> Aires écrites  : %d champs.\n', n_written);
fprintf('>> Paires ignorées (seuil manquant) : %d\n', n_skipped);
if ~isempty(notes)
    fprintf('>> Notes (%d):\n', numel(notes));
    for k = 1:numel(notes)
        fprintf('   - %s\n', notes{k});
    end
end

%% ===================== SAUVEGARDE (mise à jour en place) ================

% Réécriture du fichier .mat d'origine avec la structure DATA enrichie
save(infile, 'DATA', '-v7.3');
fprintf('>> Sauvegardé: %s\n', infile);