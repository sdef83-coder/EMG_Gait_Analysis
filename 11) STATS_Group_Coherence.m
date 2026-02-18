%% AGRÉGATION DES AIRES DE COHÉRENCE PAR GROUPE D'ÂGE — TABLES ET FIGURES
%
% OBJECTIF :
%   Ce script centralise les aires significatives de cohérence EMG-EMG
%   (CoherenceArea_*, calculées par le script 9) pour tous les participants
%   d'un projet, les organise en tables par groupe d'âge, et génère des
%   figures de barres groupées avec points individuels.
%
%   Le traitement est divisé en deux grandes parties indépendantes :
%     Partie 1 — Agrégation des données :
%       Parcourt tous les fichiers Coherence_<PID>.mat, extrait les aires
%       par condition/sous-phase/bande/paire/côté, les stocke dans ComboStore,
%       compte le nombre de cycles (N_cycle) et de participants (N_Participant)
%       par groupe/paire/condition, et exporte des CSVs + un .mat de synthèse.
%     Partie 2 — Visualisation :
%       Charge TABLES_AREA (produit en Partie 1) et génère une figure par
%       combinaison sous-phase × bande fréquentielle, avec barres groupées
%       par condition, écart-types, et points individuels par côté.
%
% ENTRÉES :
%   - Coherence_<PID>.mat : Fichiers de résultats issus des scripts 7 + 9,
%                           contenant la structure DATA avec les champs
%                           CoherenceArea_* et leurs métadonnées.
%   - ParticipantGroup.m  : Script définissant la structure Group avec les
%                           listes de participants par groupe d'âge :
%                             Group.JeunesEnfants, Group.Enfants,
%                             Group.Adolescents, Group.Adultes
%
% SORTIES :
%   Partie 1 :
%   - TABLE_AREA_<sp>_<band>_<pair>_<side>.csv : Un CSV par combinaison
%     sous-phase × bande × paire × côté, avec colonnes :
%     GroupeAge | Participant | Plat | Medium | High
%   - ALL_TABLES_AREA_STRUCT.mat : Structure TABLES_AREA complète incluant :
%       * TABLES_AREA.(pair).(band).(sp).(side) : tables MATLAB
%       * TABLES_AREA.N_cycle.(pair).(group).(cond) : total de cycles
%       * TABLES_AREA.N_Participant.(group).(pair).(cond) : nb de participants
%   Partie 2 :
%   - GroupedBars_AREA_<sp>_<band>_<groupe>.png : Une figure par combinaison
%     sous-phase × bande pour le groupe d'âge sélectionné.
%
% DÉPENDANCES :
%   - ParticipantGroup.m   : Définit la variable struct Group
%   - Fonctions locales    : resolveAge, get_quota_if_condition_has_cycles,
%                            local_side_has_cycles (en bas de ce script)
% -------------------------------------------------------------------------

%% ===================== PARTIE 1 : AGRÉGATION DES DONNÉES ================
clc; clear; close all;

%% --- Chemins ---
root_coh = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(root_coh, 'ALL');              % Fichiers Coherence_*.mat
stat_dir = fullfile(root_coh, 'STATISTIQUE_AREA'); % Dossier de sortie
csv_dir  = fullfile(stat_dir, 'CSV_CoherenceArea');
addpath(genpath('C:\Users\defsil00\Documents\Script'));

if ~exist(stat_dir,'dir'), mkdir(stat_dir); end
if ~exist(csv_dir,'dir'),  mkdir(csv_dir);  end

% Chargement des listes de participants par groupe d'âge
% → Crée la variable struct Group (JeunesEnfants, Enfants, Adolescents, Adultes)
ParticipantGroup;

% --- Paramètres ---
Conditions = {'Plat','Medium','High'};
Sides      = {'left','right','mean'};
Bandes     = {'Alpha','Beta','Gamma'};

% Paires de muscles à analyser (noms correspondant aux suffixes de champs DATA)
muscle_pairs = { ...
    'TAprox','TAdist'; ...
    'VL','RF';         ...
    'GM','SOL';        ...
    'GMED','RF';       ...
    'GMED','VL';       ...
    'RF','ST'          ...
};
% Construction des clés de paires sous la forme 'muscle1_muscle2'
pair_names = cellfun(@(a,b)[a '_' b], muscle_pairs(:,1), muscle_pairs(:,2), ...
                     'UniformOutput', false);

% Sous-phases d'intérêt par paire (détermine quels champs CoherenceArea sont extraits)
% Si une paire est absente, 'Full' est utilisé par défaut
muscle_phases = struct();
muscle_phases.TAprox_TAdist = {'LoadingResponse','Swing'};
muscle_phases.VL_RF         = {'LoadingResponse'};
muscle_phases.GM_SOL        = {'MidStance'};
muscle_phases.GMED_RF       = {'LoadingResponse'};
muscle_phases.GMED_VL       = {'LoadingResponse'};
muscle_phases.RF_ST         = {'LoadingResponse'};

% --- Détection des fichiers participants ---
mat_files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(mat_files)
    error('Aucun fichier "Coherence_*.mat" trouvé dans %s', all_dir);
end
nP = numel(mat_files); % Nombre total de participants détectés

%% --- Pré-allocation de ComboStore ---
% ComboStore est un conteneur imbriqué :
%   ComboStore.(sp).(band).(pair).(side)
% Chaque feuille contient des vecteurs de taille [nP × 1] :
%   .Participant  (cell), .GroupeAge (cell), .Plat, .Medium, .High (double)
% La pré-allocation évite la fragmentation mémoire lors de l'itération
ComboStore = struct();

for iPair = 1:numel(pair_names)
    pair = pair_names{iPair};
    if isfield(muscle_phases, pair)
        phases_for_pair = muscle_phases.(pair);
    else
        phases_for_pair = {'Full'};
    end
    for iB = 1:numel(Bandes)
        band = Bandes{iB};
        for iSP = 1:numel(phases_for_pair)
            sp = phases_for_pair{iSP};
            for iS = 1:numel(Sides)
                side = Sides{iS};
                if ~isfield(ComboStore, sp),                 ComboStore.(sp) = struct();                   end
                if ~isfield(ComboStore.(sp), band),          ComboStore.(sp).(band) = struct();             end
                if ~isfield(ComboStore.(sp).(band), pair),   ComboStore.(sp).(band).(pair) = struct();      end
                % Pré-allocation des vecteurs à NaN/cellule vide
                ComboStore.(sp).(band).(pair).(side).Participant = cell(nP,1);
                ComboStore.(sp).(band).(pair).(side).GroupeAge   = cell(nP,1);
                ComboStore.(sp).(band).(pair).(side).Plat        = nan(nP,1);
                ComboStore.(sp).(band).(pair).(side).Medium      = nan(nP,1);
                ComboStore.(sp).(band).(pair).(side).High        = nan(nP,1);
            end
        end
    end
end

% --- Initialisations pour le comptage de cycles ---
% N_cycle_total.(pair).(group).(cond) : somme des L_target_eff sur tous les
% participants du groupe (utilisée pour le panneau d'info dans les figures)
group_list    = fieldnames(Group);
N_cycle_total = struct();
counted       = struct(); % Registre anti-doublon par participant/paire/condition

% ===================== ÉTAPE 1 : Récolte des données ====================
fprintf('=== ÉTAPE 1: Harvesting AREA values + N_cycle ===\n');

for ip = 1:nP
    fn  = mat_files(ip).name;
    pid = regexprep(fn, '^Coherence_|\.mat$', ''); % Extrait l'ID du participant
    S   = load(fullfile(all_dir, fn), 'DATA');
    if ~isfield(S,'DATA'), continue; end
    DATA = S.DATA;

    grp_age = resolveAge(pid, Group); % Détermine le groupe d'âge du participant
    fprintf('Participant %d/%d: %s (%s)\n', ip, nP, pid, grp_age);

    % --- Sous-étape 1a : Comptage des cycles par paire/groupe/condition ---
    for iPair = 1:numel(pair_names)
        pair = pair_names{iPair};

        % Initialisation des nœuds si nécessaires
        if ~isfield(N_cycle_total, pair),                N_cycle_total.(pair) = struct();               end
        if ~isfield(N_cycle_total.(pair), grp_age),      N_cycle_total.(pair).(grp_age) = struct();     end
        if ~isfield(counted, pair),                      counted.(pair) = struct();                      end
        if ~isfield(counted.(pair), grp_age),            counted.(pair).(grp_age) = struct();            end

        for ic = 1:numel(Conditions)
            cond = Conditions{ic};

            if ~isfield(N_cycle_total.(pair).(grp_age), cond)
                N_cycle_total.(pair).(grp_age).(cond) = 0;
            end
            if ~isfield(counted.(pair).(grp_age), cond)
                counted.(pair).(grp_age).(cond) = struct();
            end
            if ~isfield(counted.(pair).(grp_age).(cond), pid)
                counted.(pair).(grp_age).(cond).(pid) = false;
            end

            % Chaque participant n'est compté qu'une fois par paire/cond
            if ~counted.(pair).(grp_age).(cond).(pid)
                % Lit L_target_eff uniquement si la condition a des cycles
                n_here = get_quota_if_condition_has_cycles(DATA, cond, pair);
                if ~isscalar(n_here) || ~isfinite(n_here), n_here = 0; end
                N_cycle_total.(pair).(grp_age).(cond) = ...
                    N_cycle_total.(pair).(grp_age).(cond) + double(n_here);
                counted.(pair).(grp_age).(cond).(pid) = true;
            end
        end
    end

    % --- Sous-étape 1b : Extraction des valeurs d'aire de cohérence ---
    % Seules les valeurs numériques sont stockées ici.
    % Les métadonnées Participant/GroupeAge sont remplies séparément (Étape 2)
    % pour éviter une écriture vectorielle incorrecte de cell arrays.
    for iPair = 1:numel(pair_names)
        pair  = pair_names{iPair};
        parts = strsplit(pair,'_'); m1 = parts{1}; m2 = parts{2};
        if isfield(muscle_phases, pair)
            phases_for_pair = muscle_phases.(pair);
        else
            phases_for_pair = {'Full'};
        end

        for ic = 1:numel(Conditions)
            cond = Conditions{ic};
            if ~isfield(DATA, cond), continue; end

            for iB = 1:numel(Bandes)
                band = Bandes{iB};
                for iSP = 1:numel(phases_for_pair)
                    sp = phases_for_pair{iSP};

                    % Construction du nom du champ CoherenceArea selon la sous-phase
                    if strcmpi(sp,'Full')
                        area_field = sprintf('CoherenceArea_%s_%s_%s', band, m1, m2);
                    else
                        area_field = sprintf('CoherenceArea_%s_%s_%s_%s', sp, band, m1, m2);
                    end

                    % Lecture des valeurs gauche, droite et calcul de la moyenne
                    vL = NaN; vR = NaN;
                    if isfield(DATA.(cond),'left')  && isfield(DATA.(cond).left,  area_field)
                        vL = DATA.(cond).left.(area_field);
                    end
                    if isfield(DATA.(cond),'right') && isfield(DATA.(cond).right, area_field)
                        vR = DATA.(cond).right.(area_field);
                    end
                    vM = mean([vL, vR], 'omitnan'); % Moyenne bilatérale

                    % Stockage dans le vecteur pré-alloué à l'indice ip du participant
                    switch cond
                        case 'Plat'
                            ComboStore.(sp).(band).(pair).left.Plat(ip)   = vL;
                            ComboStore.(sp).(band).(pair).right.Plat(ip)  = vR;
                            ComboStore.(sp).(band).(pair).mean.Plat(ip)   = vM;
                        case 'Medium'
                            ComboStore.(sp).(band).(pair).left.Medium(ip)  = vL;
                            ComboStore.(sp).(band).(pair).right.Medium(ip) = vR;
                            ComboStore.(sp).(band).(pair).mean.Medium(ip)  = vM;
                        case 'High'
                            ComboStore.(sp).(band).(pair).left.High(ip)   = vL;
                            ComboStore.(sp).(band).(pair).right.High(ip)  = vR;
                            ComboStore.(sp).(band).(pair).mean.High(ip)   = vM;
                    end
                end
            end
        end
    end
end

% ===================== ÉTAPE 2 : Remplissage des métadonnées ============
% Séparé de l'Étape 1 pour éviter que l'écriture dans des cell arrays
% imbriqués n'écrase des entrées voisines (chaque indice ip est isolé)
fprintf('=== ÉTAPE 2: Remplissage des métadonnées Participant/GroupeAge ===\n');

for ip = 1:nP
    fn      = mat_files(ip).name;
    pid     = regexprep(fn, '^Coherence_|\.mat$', '');
    grp_age = resolveAge(pid, Group);
    fprintf('Métadonnées participant %d/%d: %s (%s)\n', ip, nP, pid, grp_age);

    % Parcours exhaustif de ComboStore pour remplir uniquement l'index ip
    sp_fields = fieldnames(ComboStore);
    for isf = 1:numel(sp_fields)
        sp = sp_fields{isf};
        band_fields = fieldnames(ComboStore.(sp));
        for ibf = 1:numel(band_fields)
            band = band_fields{ibf};
            pair_fields = fieldnames(ComboStore.(sp).(band));
            for ipf = 1:numel(pair_fields)
                pair = pair_fields{ipf};
                sides_fields = fieldnames(ComboStore.(sp).(band).(pair));
                for isd = 1:numel(sides_fields)
                    side = sides_fields{isd};
                    ComboStore.(sp).(band).(pair).(side).Participant{ip} = pid;
                    ComboStore.(sp).(band).(pair).(side).GroupeAge{ip}   = grp_age;
                end
            end
        end
    end
end

% ===================== ÉTAPE 3 : Construction des tables et CSV =========
% Pour chaque combinaison (sp, band, pair, side), construit une table MATLAB
% et exporte un CSV dans csv_dir
TABLES_AREA = struct();
fprintf('=== ÉTAPE 3: Writing CSVs and assembling TABLES_AREA ===\n');

sp_fields = fieldnames(ComboStore);
for isf = 1:numel(sp_fields)
    sp = sp_fields{isf};
    band_fields = fieldnames(ComboStore.(sp));
    for ibf = 1:numel(band_fields)
        band = band_fields{ibf};
        pair_fields = fieldnames(ComboStore.(sp).(band));
        for ipf = 1:numel(pair_fields)
            pair = pair_fields{ipf};
            sides_fields = fieldnames(ComboStore.(sp).(band).(pair));

            % Initialisation de la branche dans TABLES_AREA
            if ~isfield(TABLES_AREA, pair),              TABLES_AREA.(pair) = struct();              end
            if ~isfield(TABLES_AREA.(pair), band),       TABLES_AREA.(pair).(band) = struct();       end
            if ~isfield(TABLES_AREA.(pair).(band), sp),  TABLES_AREA.(pair).(band).(sp) = struct();  end

            % Un CSV par côté (left, right, mean)
            for isd = 1:numel(sides_fields)
                side  = sides_fields{isd};
                store = ComboStore.(sp).(band).(pair).(side);

                T = table( ...
                    store.GroupeAge, store.Participant, ...
                    store.Plat, store.Medium, store.High, ...
                    'VariableNames', {'GroupeAge','Participant','Plat','Medium','High'});

                TABLES_AREA.(pair).(band).(sp).(side) = T;

                csv_name = sprintf('TABLE_AREA_%s_%s_%s_%s.csv', sp, band, pair, side);
                writetable(T, fullfile(csv_dir, csv_name));
            end
        end
    end
end

% --- Vérification rapide sur une combinaison connue ---
fprintf('=== ÉTAPE 4: Vérification rapide ===\n');
if isfield(TABLES_AREA, 'TAprox_TAdist') && ...
   isfield(TABLES_AREA.TAprox_TAdist, 'Alpha') && ...
   isfield(TABLES_AREA.TAprox_TAdist.Alpha, 'LoadingResponse') && ...
   isfield(TABLES_AREA.TAprox_TAdist.Alpha.LoadingResponse, 'mean')

    T_test = TABLES_AREA.TAprox_TAdist.Alpha.LoadingResponse.mean;
    fprintf('Test sur TAprox_TAdist/Alpha/LoadingResponse/mean:\n');
    for i = 1:min(5, height(T_test))
        fprintf('  Participant %d: %s (%s) - Plat: %.3f\n', i, ...
            T_test.Participant{i}, T_test.GroupeAge{i}, T_test.Plat(i));
    end
else
    fprintf('Structure de test non trouvée - vérifiez vos données\n');
end

% ===================== ÉTAPE 4 (suite) : Injection de N_cycle ===========
% Ajoute le total de cycles par paire/groupe/condition dans TABLES_AREA
% pour affichage dans le panneau d'info des figures
TABLES_AREA.N_cycle = struct();

pairs = fieldnames(N_cycle_total);
for i = 1:numel(pairs)
    pair   = pairs{i};
    groups = fieldnames(N_cycle_total.(pair));
    for g = 1:numel(groups)
        grp = groups{g};
        if ~isfield(TABLES_AREA.N_cycle, pair),        TABLES_AREA.N_cycle.(pair) = struct();        end
        if ~isfield(TABLES_AREA.N_cycle.(pair), grp),  TABLES_AREA.N_cycle.(pair).(grp) = struct();  end

        for ic = 1:numel(Conditions)
            cond = Conditions{ic};
            val = 0;
            if isfield(N_cycle_total.(pair).(grp), cond)
                val = N_cycle_total.(pair).(grp).(cond);
            end
            TABLES_AREA.N_cycle.(pair).(grp).(cond) = val;
        end
    end
end

% ===================== ÉTAPE 5 : Comptage des participants ==============
% Pour chaque combinaison groupe/paire/condition, compte le nombre de
% participants ayant au moins une valeur non-NaN (contribution réelle)
fprintf('=== Computing N_Participant by group/pair/condition ===\n');

if ~isfield(TABLES_AREA,'N_Participant'), TABLES_AREA.N_Participant = struct(); end

all_pairs = fieldnames(TABLES_AREA);
all_pairs = setdiff(all_pairs, {'N_cycle','N_Participant'}); % Exclut les branches de synthèse

conds  = {'Plat','Medium','High'};
groups = fieldnames(Group);

for ipair = 1:numel(all_pairs)
    pair = all_pairs{ipair};

    % contrib.(group).(pid) = [logical x 3] : a une valeur valide dans [Plat, Medium, High]
    contrib = struct();
    for ig = 1:numel(groups), contrib.(groups{ig}) = struct(); end

    if ~isstruct(TABLES_AREA.(pair)), continue; end
    bands = fieldnames(TABLES_AREA.(pair));
    for ib = 1:numel(bands)
        band = bands{ib};
        if ~isstruct(TABLES_AREA.(pair).(band)), continue; end
        subphases = fieldnames(TABLES_AREA.(pair).(band));
        for isp = 1:numel(subphases)
            sp = subphases{isp};
            if ~isfield(TABLES_AREA.(pair).(band).(sp), 'mean'), continue; end
            T = TABLES_AREA.(pair).(band).(sp).mean;
            if isempty(T), continue; end

            % Pour chaque ligne (participant), enregistre s'il contribue
            % (valeur non-NaN) à chaque condition
            for irow = 1:height(T)
                gr  = string(T.GroupeAge{irow});
                pid = string(T.Participant{irow});
                if strlength(gr)==0 || strlength(pid)==0, continue; end
                gr = char(gr); pid = char(pid);

                if ~isfield(contrib.(gr), pid)
                    contrib.(gr).(pid) = false(1, numel(conds));
                end
                vals = [T.Plat(irow), T.Medium(irow), T.High(irow)];
                contrib.(gr).(pid) = contrib.(gr).(pid) | ~isnan(vals);
            end
        end
    end

    % Compte par groupe et condition en n'incluant que les participants
    % ayant contribué à AU MOINS une condition (exclut les lignes vides)
    for ig = 1:numel(groups)
        gr   = groups{ig};
        pids = fieldnames(contrib.(gr));
        if isempty(pids)
            for ic = 1:numel(conds)
                cond = conds{ic};
                if ~isfield(TABLES_AREA.N_Participant, gr),         TABLES_AREA.N_Participant.(gr) = struct();          end
                if ~isfield(TABLES_AREA.N_Participant.(gr), pair),  TABLES_AREA.N_Participant.(gr).(pair) = struct();   end
                TABLES_AREA.N_Participant.(gr).(pair).(cond) = 0;
            end
            continue;
        end

        % Matrice [nParticipants × 3 conditions] de contributions
        M = false(numel(pids), numel(conds));
        for k = 1:numel(pids), M(k,:) = contrib.(gr).(pids{k}); end
        M = M(any(M,2), :); % Garde uniquement les participants actifs

        for ic = 1:numel(conds)
            cond = conds{ic};
            n_here = sum(M(:,ic));
            if ~isfield(TABLES_AREA.N_Participant, gr),        TABLES_AREA.N_Participant.(gr) = struct();         end
            if ~isfield(TABLES_AREA.N_Participant.(gr), pair), TABLES_AREA.N_Participant.(gr).(pair) = struct();  end
            TABLES_AREA.N_Participant.(gr).(pair).(cond) = n_here;
        end
    end
end

% ===================== SAUVEGARDE DU .mat MAÎTRE ========================
mat_out = fullfile(stat_dir, 'ALL_TABLES_AREA_STRUCT.mat');
save(mat_out, 'TABLES_AREA', '-v7.3');
fprintf('DONE: %s\n', mat_out);

%% ===================== PARTIE 2 : VISUALISATION =========================
clc; close all;

% --- Chemins ---
root_coh = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(root_coh, 'ALL');
stat_dir = fullfile(root_coh, 'STATISTIQUE_AREA');
csv_dir  = fullfile(stat_dir, 'CSV_CoherenceArea');
addpath(genpath('C:\Users\defsil00\Documents\Script'));

% --- Paramètres de visualisation ---
Conditions = {'Plat','Medium','High'};
Bandes     = {'Alpha','Beta','Gamma'};

% Plages fréquentielles pour les étiquettes de figures
band_ranges = struct('Alpha','8-12 Hz', 'Beta','13-30 Hz', 'Gamma','31-60 Hz');

% Esthétique des figures
col_left   = [0.00 0.60 0.00]; % Vert — points côté gauche
col_right  = [0.50 0.00 0.50]; % Violet — points côté droit
alpha_pts  = 0.40;              % Transparence des points individuels
% Niveaux de gris pour les barres des 3 conditions (Plat → plus clair)
gray_levels = [0.85 0.85 0.85; 0.70 0.70 0.70; 0.55 0.55 0.55];

% Dossier de sortie des figures
fig_dir = fullfile(stat_dir, 'Fig_CoherenceArea_GroupedBars');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Chargement de TABLES_AREA si non présent en workspace
if ~exist('TABLES_AREA','var')
    S = load(fullfile(stat_dir,'ALL_TABLES_AREA_STRUCT.mat'));
    TABLES_AREA = S.TABLES_AREA;
end

% Chargement des groupes d'âge
ParticipantGroup;

% Groupe d'âge à visualiser (adapter selon besoin)
groupe_age = 'Adolescents'; % 'JeunesEnfants' | 'Enfants' | 'Adolescents' | 'Adultes'

% Extraction des paires disponibles (hors branches de synthèse N_cycle, N_Participant)
all_pairs = setdiff(fieldnames(TABLES_AREA), {'N_cycle','N_Participant'});

% Inventaire de toutes les sous-phases effectivement présentes dans TABLES_AREA
all_subphases = {};
for ip = 1:numel(all_pairs)
    p = all_pairs{ip};
    bnames = fieldnames(TABLES_AREA.(p));
    for ib = 1:numel(bnames)
        b = bnames{ib};
        spnames = fieldnames(TABLES_AREA.(p).(b));
        all_subphases = [all_subphases, spnames']; 
    end
end
all_subphases = unique(all_subphases);

fprintf('=== Grouped bar plots for COHERENCE AREA (groupe: %s) ===\n', groupe_age);

% --- Boucle de génération des figures ---
for ib = 1:numel(Bandes)
    band = Bandes{ib};

    for isp = 1:numel(all_subphases)
        sp = all_subphases{isp};

        % Sélection des paires ayant des données pour cette combinaison band/sp
        pair_list = {};
        for ip = 1:numel(all_pairs)
            p = all_pairs{ip};
            if isfield(TABLES_AREA.(p), band) && isfield(TABLES_AREA.(p).(band), sp) && ...
               all(isfield(TABLES_AREA.(p).(band).(sp), {'left','right','mean'}))
                pair_list{end+1} = p;
            end
        end
        if isempty(pair_list)
            fprintf('  [!] Aucune paire pour subphase %s / bande %s\n', sp, band);
            continue;
        end

        % Calcul des moyennes et écarts-types par condition et paire
        % (filtrés sur le groupe d'âge sélectionné)
        M  = nan(numel(Conditions), numel(pair_list)); % Moyennes
        SD = nan(numel(Conditions), numel(pair_list)); % Écarts-types

        % Stockage des tables filtrées pour les scatters individuels
        Tm_cell = cell(numel(pair_list),1);
        Tl_cell = cell(numel(pair_list),1);
        Tr_cell = cell(numel(pair_list),1);

        for ip = 1:numel(pair_list)
            p  = pair_list{ip};
            Tm = TABLES_AREA.(p).(band).(sp).mean;
            Tl = TABLES_AREA.(p).(band).(sp).left;
            Tr = TABLES_AREA.(p).(band).(sp).right;

            % Filtrage par groupe d'âge
            Tm = Tm(ismember(Tm.GroupeAge, groupe_age), :);
            Tl = Tl(ismember(Tl.GroupeAge, groupe_age), :);
            Tr = Tr(ismember(Tr.GroupeAge, groupe_age), :);
            Tm_cell{ip} = Tm; Tl_cell{ip} = Tl; Tr_cell{ip} = Tr;

            for ic = 1:numel(Conditions)
                c = Conditions{ic};
                vals       = Tm.(c);
                M(ic, ip)  = mean(vals, 'omitnan');
                SD(ic, ip) = std(vals,  'omitnan');
            end
        end

        % --- Tracé de la figure ---
        x  = 1:numel(pair_list);
        f  = figure('Position',[100 120 1600 620]); hold on;

        % Barres groupées (conditions = groupes, paires = positions sur x)
        bh = bar(x, M', 'grouped'); % M' = [nPairs × nConditions]
        for ic = 1:numel(Conditions)
            bh(ic).FaceColor = gray_levels(ic,:);
            bh(ic).EdgeColor = [0 0 0];
        end
        drawnow;

        % Récupération des positions exactes du centre de chaque barre
        xC = arrayfun(@(b)b.XEndPoints, bh, 'UniformOutput', false);

        % Superposition des points de moyenne et des barres d'erreur (écart-type)
        for ic = 1:numel(Conditions)
            xc = xC{ic};
            for ip = 1:numel(pair_list)
                mu = M(ic,ip); sd = SD(ic,ip);
                if ~isnan(mu)
                    plot(xc(ip), mu, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
                    errorbar(xc(ip), mu, sd, 'k', 'LineWidth',1.2, 'CapSize',8);
                end
            end
        end

        % Superposition des points individuels par côté avec jitter horizontal
        % (améliore la lisibilité en évitant la superposition des points)
        jitter = 0.06;
        for ic = 1:numel(Conditions)
            xc = xC{ic};
            for ip = 1:numel(pair_list)
                % Points côté gauche (vert)
                if ~isempty(Tl_cell{ip})
                    vL = Tl_cell{ip}.(Conditions{ic});
                    maskL = ~isnan(vL);
                    if any(maskL)
                        xL = xc(ip) - jitter + (rand(sum(maskL),1)*2*jitter);
                        sL = scatter(xL, vL(maskL), 28, 'o', ...
                                     'MarkerFaceColor',col_left, 'MarkerEdgeColor',col_left);
                        sL.MarkerFaceAlpha = alpha_pts; sL.MarkerEdgeAlpha = alpha_pts;
                    end
                end
                % Points côté droit (violet)
                if ~isempty(Tr_cell{ip})
                    vR = Tr_cell{ip}.(Conditions{ic});
                    maskR = ~isnan(vR);
                    if any(maskR)
                        xR = xc(ip) - jitter + (rand(sum(maskR),1)*2*jitter);
                        sR = scatter(xR, vR(maskR), 28, 'o', ...
                                     'MarkerFaceColor',col_right, 'MarkerEdgeColor',col_right);
                        sR.MarkerFaceAlpha = alpha_pts; sR.MarkerEdgeAlpha = alpha_pts;
                    end
                end
            end
        end

        % Formatage des axes
        set(gca, 'XTick', x, 'XTickLabel', strrep(pair_list,'_','-'), ...
                 'TickLabelInterpreter','none');
        xlabel('Paires musculaires', 'FontSize',12, 'FontWeight','bold');
        ylabel('Coherence area (a.u.)', 'FontSize',12, 'FontWeight','bold');

        % Légende
        hLeft  = scatter(nan, nan, 36, 'o', 'MarkerFaceColor',col_left,  'MarkerEdgeColor',col_left);
        hRight = scatter(nan, nan, 36, 'o', 'MarkerFaceColor',col_right, 'MarkerEdgeColor',col_right);
        legend([bh(1), bh(2), bh(3), hLeft, hRight], ...
               {'Plat','Medium','High','left','right'}, 'Location','northwest');

        % --- Panneau d'information latéral (N participants + N cycles) ---
        % Colonnes : PAIRE | N_Plat N_Medium N_High | C_Plat C_Medium C_High
        wPair = 14; wNs = 5; wCs = 6;
        fmtHeader = sprintf('%%-%ds  %%%ds %%%ds %%%ds   %%%ds %%%ds %%%ds', wPair, wNs, wNs, wNs, wCs, wCs, wCs);
        fmtRow    = sprintf('%%-%ds  %%%dd %%%dd %%%dd   %%%dd %%%dd %%%dd', wPair, wNs, wNs, wNs, wCs, wCs, wCs);

        info_lines = cell(1, numel(pair_list)+1);
        info_lines{1} = sprintf(fmtHeader, 'PAIRE', 'N_P', 'N_M', 'N_H', 'C_P', 'C_M', 'C_H');

        for ipx = 1:numel(pair_list)
            pp      = strrep(string(pair_list{ipx}), '_', '-');
            Tm_pair = Tm_cell{ipx};

            % N participants avec valeur non-NaN par condition (table filtrée au groupe)
            Np = [0 0 0];
            if ~isempty(Tm_pair)
                for ic = 1:numel(Conditions)
                    c = Conditions{ic};
                    if ismember(c, Tm_pair.Properties.VariableNames)
                        Np(ic) = sum(~isnan(Tm_pair.(c)));
                    end
                end
            end

            % C total de cycles (depuis TABLES_AREA.N_cycle)
            Cp = [0 0 0];
            for ic = 1:numel(Conditions)
                cond = Conditions{ic};
                if isfield(TABLES_AREA,'N_cycle') && ...
                   isfield(TABLES_AREA.N_cycle, pair_list{ipx}) && ...
                   isfield(TABLES_AREA.N_cycle.(pair_list{ipx}), groupe_age) && ...
                   isfield(TABLES_AREA.N_cycle.(pair_list{ipx}).(groupe_age), cond)
                    Cp(ic) = TABLES_AREA.N_cycle.(pair_list{ipx}).(groupe_age).(cond);
                end
            end
            info_lines{ipx+1} = sprintf(fmtRow, char(pp), Np(1), Np(2), Np(3), Cp(1), Cp(2), Cp(3));
        end

        % Positionnement de l'axe principal pour laisser de la place au panneau
        ax = gca;
        ax.Position = [0.08 0.18 0.56 0.72];

        % Titre
        t = title(sprintf('%s (%s) - %s - %s', band, band_ranges.(band), sp, groupe_age), ...
                  'FontSize',14,'FontWeight','bold');
        t.Units = 'normalized'; t.Position(2) = 1.03;

        % Annotation textbox (police Consolas pour l'alignement en colonnes)
        annotation(f, 'textbox', [0.68 0.18 0.32 0.72], ...
                   'String', strjoin(info_lines, newline), ...
                   'Interpreter','none', ...
                   'FontName','Consolas', 'FontSize',10, ...
                   'HorizontalAlignment','left', ...
                   'VerticalAlignment','top', ...
                   'EdgeColor',[0.85 0.85 0.85], 'BackgroundColor',[1 1 1], ...
                   'FitBoxToText','off', 'Margin', 6);

        % Export PNG haute résolution
        figname = sprintf('GroupedBars_AREA_%s_%s_%s.png', sp, band, groupe_age);
        exportgraphics(f, fullfile(fig_dir, figname), 'Resolution', 300);
        fprintf('  OK FIG: %s / %s -> %s\n', sp, band, figname);
        close(f);
    end
end

fprintf('All done.\n');

%% ===================== FONCTIONS LOCALES ================================

function grp = resolveAge(pid, Group)
% RESOLVEAGE  Retourne le nom du groupe d'âge d'un participant.
%
%   Entrée  : pid   — identifiant du participant (string)
%             Group — struct avec champs JeunesEnfants, Enfants,
%                     Adolescents (optionnel), Adultes (listes de PIDs)
%   Sortie  : grp   — nom du groupe ('JeunesEnfants','Enfants',
%                     'Adolescents','Adultes', ou 'Inconnu')

    if ismember(pid, Group.JeunesEnfants)
        grp = 'JeunesEnfants';
    elseif ismember(pid, Group.Enfants)
        grp = 'Enfants';
    elseif isfield(Group,'Adolescents') && ismember(pid, Group.Adolescents)
        grp = 'Adolescents';
    elseif ismember(pid, Group.Adultes)
        grp = 'Adultes';
    else
        grp = 'Inconnu';
    end
end

function n = get_quota_if_condition_has_cycles(DATA, cond, pair)
% GET_QUOTA_IF_CONDITION_HAS_CYCLES
%   Retourne L_target_eff_<pair> pour la condition donnée, mais uniquement
%   si au moins un côté (gauche ou droit) possède des cycles retenus.
%   Retourne 0 si la condition est vide ou si le quota est introuvable.
%
%   Entrées : DATA  — structure DATA du participant
%             cond  — nom de la condition ('Plat', 'Medium', 'High')
%             pair  — clé de paire ('m1_m2')
%   Sortie  : n     — nombre de cycles (scalaire double), ou 0

    n = 0;
    if ~isstruct(DATA) || ~isfield(DATA, cond), return; end
    parts = strsplit(pair,'_');
    if numel(parts) ~= 2, return; end
    m1 = parts{1}; m2 = parts{2};

    % Vérifie la présence de cycles sur au moins un côté
    has_cycles = local_side_has_cycles(DATA.(cond), 'left',  m1, m2) || ...
                 local_side_has_cycles(DATA.(cond), 'right', m1, m2);
    if ~has_cycles, return; end

    % Lecture du quota dans les métadonnées de la condition
    if isfield(DATA.(cond), 'meta') && isstruct(DATA.(cond).meta)
        meta = DATA.(cond).meta;
        f1 = ['L_target_eff_' pair];
        f2 = ['L_target_'     pair]; % Fallback pour anciens fichiers
        if isfield(meta, f1) && isnumeric(meta.(f1)) && isfinite(meta.(f1))
            n = double(meta.(f1)); return;
        elseif isfield(meta, f2) && isnumeric(meta.(f2)) && isfinite(meta.(f2))
            n = double(meta.(f2)); return;
        end
    end
    n = 0; % Quota introuvable dans les métadonnées
end

function tf = local_side_has_cycles(node_cond, side, m1, m2)
% LOCAL_SIDE_HAS_CYCLES
%   Retourne true si le côté donné possède au moins 1 cycle retenu
%   pour la paire m1-m2, en vérifiant les champs Ncycle_<m1>_<m2>
%   puis L_<m1>_<m2> (fallback pour compatibilité avec anciens fichiers).
%
%   Entrées : node_cond — DATA.(cond)
%             side      — 'left' ou 'right'
%             m1, m2    — noms des muscles de la paire
%   Sortie  : tf        — booléen

    tf = false;
    if ~isfield(node_cond, side), return; end
    S = node_cond.(side);

    % Vérification via Ncycle_<m1>_<m2> (champ principal)
    fnN = ['Ncycle_' m1 '_' m2];
    if isfield(S, fnN) && isnumeric(S.(fnN)) && isfinite(S.(fnN)) && S.(fnN) > 0
        tf = true; return;
    end
    % Fallback : vérification via L_<m1>_<m2> (anciens fichiers)
    fnL = ['L_' m1 '_' m2];
    if isfield(S, fnL) && isnumeric(S.(fnL)) && isfinite(S.(fnL)) && S.(fnL) > 0
        tf = true; return;
    end
end