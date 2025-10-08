%% ========================================================================
%  COHERENCE AREA - SOLUTION 2: SÉPARATION COMPLÈTE DES REMPLISSAGES
% ========================================================================

%=========================== 1) ORGANIZATION =============================
clc; clear; close all;

% --- Folders (EDIT IF NEEDED) -------------------------------------------
root_coh = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(root_coh, 'ALL');                 % contains Coherence_*.mat
stat_dir = fullfile(root_coh, 'STATISTIQUE_AREA');    % new output root for AREA
csv_dir  = fullfile(stat_dir, 'CSV_CoherenceArea');
addpath(genpath('C:\Users\defsil00\Documents\Script'));

if ~exist(stat_dir,'dir'), mkdir(stat_dir); end
if ~exist(csv_dir,'dir'), mkdir(csv_dir); end

% --- Add utils path (EDIT IF NEEDED) ------------------------------------
addpath(genpath('C:\Users\defsil00\Documents\Script'));

% --- Load age groups definition -----------------------------------------
% Must define struct "Group" with fields: JeunesEnfants, Enfants, Adolescents (optional), Adultes
ParticipantGroup;

% --- Parameters ----------------------------------------------------------
Conditions = {'Plat','Medium','High'};
Sides      = {'left','right','mean'};
Bandes     = {'Alpha','Beta','Gamma'};

% Muscle pairs (names must match your DATA field suffixes)
muscle_pairs = { ...
    'TAprox','TAdist'; ...
    'VL','RF'; ...
    'GM','SOL'; ...
    'GMED','RF'; ...
    'GMED','VL'; ...
    'RF','ST' ...
};
pair_names = cellfun(@(a,b)[a '_' b], muscle_pairs(:,1), muscle_pairs(:,2), 'UniformOutput', false);

% Which subphases to include per pair (mirrors your plotting logic)
muscle_phases = struct();
muscle_phases.TAprox_TAdist = {'LoadingResponse','Swing'};
muscle_phases.VL_RF         = {'LoadingResponse'};
muscle_phases.GM_SOL        = {'MidStance'};
muscle_phases.GMED_RF       = {'LoadingResponse'};
muscle_phases.GMED_VL       = {'LoadingResponse'};
muscle_phases.RF_ST         = {'LoadingResponse'};

% --- List participant files ---------------------------------------------
mat_files = dir(fullfile(all_dir, 'Coherence_*.mat'));
if isempty(mat_files)
    error('No "Coherence_*.mat" files found in %s', all_dir);
end
nP = numel(mat_files);

% --- Nested container: ComboStore.(sp).(band).(pair).(side) -------------
ComboStore = struct(); % will hold fields: Participant (nP x 1), GroupeAge (nP x 1), Plat/Medium/High (nP x 1)

% Pre-build all (sp, band, pair, side) slots with preallocation -----------
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
                if ~isfield(ComboStore, sp), ComboStore.(sp) = struct(); end
                if ~isfield(ComboStore.(sp), band), ComboStore.(sp).(band) = struct(); end
                if ~isfield(ComboStore.(sp).(band), pair), ComboStore.(sp).(band).(pair) = struct(); end
                % Preallocate vectors
                ComboStore.(sp).(band).(pair).(side).Participant = cell(nP,1);
                ComboStore.(sp).(band).(pair).(side).GroupeAge  = cell(nP,1);
                ComboStore.(sp).(band).(pair).(side).Plat       = nan(nP,1);
                ComboStore.(sp).(band).(pair).(side).Medium     = nan(nP,1);
                ComboStore.(sp).(band).(pair).(side).High       = nan(nP,1);
            end
        end
    end
end

% --- Totaux de cycles par paire/groupe/condition -------------------------
group_list = fieldnames(Group); % ex: {'JeunesEnfants','Enfants','Adolescents','Adultes'}

% N_cycle_total.(pair).(group).(cond) = somme des L_target_eff_<pair> lus dans META de chaque condition
N_cycle_total = struct();

% Registre "déjà compté" (évite de recompter le même participant)
% counted.(pair).(group).(cond).(pid) = true
counted = struct();

fprintf('=== ÉTAPE 1: Harvesting AREA values + N_cycle ===\n');

% ✅ BOUCLE PRINCIPALE: SEULEMENT RÉCUPÉRATION DES DONNÉES
for ip = 1:nP
    fn = mat_files(ip).name;            % Coherence_CTL_xx.mat
    pid = regexprep(fn, '^Coherence_|\.mat$', '');
    S = load(fullfile(all_dir, fn), 'DATA');
    if ~isfield(S,'DATA'), continue; end
    DATA = S.DATA;

    grp_age = resolveAge(pid, Group);
    fprintf('Participant %d/%d: %s (%s)\n', ip, nP, pid, grp_age);

    % === 1) COMPTER LES CYCLES (INCHANGÉ) ===================================
    for iPair = 1:numel(pair_names)
        pair = pair_names{iPair};

        % Init si besoin
        if ~isfield(N_cycle_total, pair), N_cycle_total.(pair) = struct(); end
        if ~isfield(N_cycle_total.(pair), grp_age), N_cycle_total.(pair).(grp_age) = struct(); end
        if ~isfield(counted, pair), counted.(pair) = struct(); end
        if ~isfield(counted.(pair), grp_age), counted.(pair).(grp_age) = struct(); end

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

            if ~counted.(pair).(grp_age).(cond).(pid)
                n_here = get_quota_if_condition_has_cycles(DATA, cond, pair);  % <-- NOUVEAU
                if ~isscalar(n_here) || ~isfinite(n_here), n_here = 0; end
                N_cycle_total.(pair).(grp_age).(cond) = N_cycle_total.(pair).(grp_age).(cond) + double(n_here);
                counted.(pair).(grp_age).(cond).(pid) = true;
            end
        end
    end

    % === 2) RÉCUPÉRER LES AIRES - SANS REMPLIR Participant/GroupeAge ======
    for iPair = 1:numel(pair_names)
        pair = pair_names{iPair};
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

                    % construire area_field
                    if strcmpi(sp,'Full')
                        area_field = sprintf('CoherenceArea_%s_%s_%s', band, m1, m2);
                    else
                        area_field = sprintf('CoherenceArea_%s_%s_%s_%s', sp, band, m1, m2);
                    end

                    % Fetch left/right/mean
                    vL = NaN; vR = NaN;
                    if isfield(DATA.(cond),'left')  && isfield(DATA.(cond).left,  area_field)
                        vL = DATA.(cond).left.(area_field);
                    end
                    if isfield(DATA.(cond),'right') && isfield(DATA.(cond).right, area_field)
                        vR = DATA.(cond).right.(area_field);
                    end
                    vM = mean([vL, vR], 'omitnan');

                    % Store SEULEMENT LES VALEURS NUMÉRIQUES
                    switch cond
                        case 'Plat'
                            ComboStore.(sp).(band).(pair).left.Plat(ip)  = vL;
                            ComboStore.(sp).(band).(pair).right.Plat(ip) = vR;
                            ComboStore.(sp).(band).(pair).mean.Plat(ip)  = vM;
                        case 'Medium'
                            ComboStore.(sp).(band).(pair).left.Medium(ip)  = vL;
                            ComboStore.(sp).(band).(pair).right.Medium(ip) = vR;
                            ComboStore.(sp).(band).(pair).mean.Medium(ip)  = vM;
                        case 'High'
                            ComboStore.(sp).(band).(pair).left.High(ip)  = vL;
                            ComboStore.(sp).(band).(pair).right.High(ip) = vR;
                            ComboStore.(sp).(band).(pair).mean.High(ip)  = vM;
                    end
                end
            end
        end
    end
end

% ✅ ÉTAPE 2: REMPLISSAGE SÉPARÉ DES MÉTADONNÉES Participant/GroupeAge
fprintf('=== ÉTAPE 2: Remplissage des métadonnées Participant/GroupeAge ===\n');

for ip = 1:nP
    fn = mat_files(ip).name;
    pid = regexprep(fn, '^Coherence_|\.mat$', '');
    grp_age = resolveAge(pid, Group);
    
    fprintf('Métadonnées participant %d/%d: %s (%s)\n', ip, nP, pid, grp_age);
    
    % Parcourir TOUTES les combinaisons et remplir SEULEMENT l'index ip
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
                    % ✅ REMPLIR SEULEMENT L'INDEX ip, PAS TOUS LES INDICES
                    ComboStore.(sp).(band).(pair).(side).Participant{ip} = pid;
                    ComboStore.(sp).(band).(pair).(side).GroupeAge{ip}   = grp_age;
                end
            end
        end
    end
end

%=================== 3) Build TABLES_AREA + CSVs (INCHANGÉ) ==============
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

            % Create branch in TABLES_AREA
            if ~isfield(TABLES_AREA, pair), TABLES_AREA.(pair) = struct(); end
            if ~isfield(TABLES_AREA.(pair), band), TABLES_AREA.(pair).(band) = struct(); end
            if ~isfield(TABLES_AREA.(pair).(band), sp), TABLES_AREA.(pair).(band).(sp) = struct(); end

            % Emit a CSV per side
            for isd = 1:numel(sides_fields)
                side = sides_fields{isd};
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

% ✅ TEST DE VÉRIFICATION
fprintf('=== ÉTAPE 4: Vérification rapide ===\n');
if isfield(TABLES_AREA, 'TAprox_TAdist') && ...
   isfield(TABLES_AREA.TAprox_TAdist, 'Alpha') && ...
   isfield(TABLES_AREA.TAprox_TAdist.Alpha, 'LoadingResponse') && ...
   isfield(TABLES_AREA.TAprox_TAdist.Alpha.LoadingResponse, 'mean')
    
    T_test = TABLES_AREA.TAprox_TAdist.Alpha.LoadingResponse.mean;
    fprintf('Test sur TAprox_TAdist/Alpha/LoadingResponse/mean:\n');
    for i = 1:min(5, height(T_test))
        pid = T_test.Participant{i};
        grp = T_test.GroupeAge{i};
        plat = T_test.Plat(i);
        fprintf('  Participant %d: %s (%s) - Plat: %.3f\n', i, pid, grp, plat);
    end
else
    fprintf('Structure de test non trouvée - vérifiez vos données\n');
end

%=================== 4) Injecter N_cycle (par condition) =================
TABLES_AREA.N_cycle = struct();

pairs = fieldnames(N_cycle_total);
for i = 1:numel(pairs)
    pair = pairs{i};
    groups = fieldnames(N_cycle_total.(pair));
    for g = 1:numel(groups)
        grp = groups{g};
        if ~isfield(TABLES_AREA.N_cycle, pair), TABLES_AREA.N_cycle.(pair) = struct(); end
        if ~isfield(TABLES_AREA.N_cycle.(pair), grp), TABLES_AREA.N_cycle.(pair).(grp) = struct(); end

        % Remplir par condition (sans copie forcée)
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

%============= 5) NOMBRE DE PARTICIPANTS PAR GROUPE/PAIRE/CONDITION ======
fprintf('=== Computing N_Participant by group/pair/condition ===\n');

if ~isfield(TABLES_AREA,'N_Participant'), TABLES_AREA.N_Participant = struct(); end

all_pairs = fieldnames(TABLES_AREA);
all_pairs = setdiff(all_pairs, {'N_cycle','N_Participant'});  % ignorer branches de synthèse

conds = {'Plat','Medium','High'};
groups = fieldnames(Group);  % mêmes groupes que ParticipantGroup.m

for ipair = 1:numel(all_pairs)
    pair = all_pairs{ipair};

    % Registre des contributions: contrib.(group).(pid) = [c1 c2 c3]
    contrib = struct();
    for ig = 1:numel(groups)
        contrib.(groups{ig}) = struct();
    end

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

            for irow = 1:height(T)
                gr  = string(T.GroupeAge{irow});
                pid = string(T.Participant{irow});
                if strlength(gr)==0 || strlength(pid)==0, continue; end
                gr  = char(gr); pid = char(pid);

                if ~isfield(contrib.(gr), pid)
                    contrib.(gr).(pid) = false(1, numel(conds));
                end

                vals = [T.Plat(irow), T.Medium(irow), T.High(irow)];
                mask_valid = ~isnan(vals);
                contrib.(gr).(pid) = contrib.(gr).(pid) | mask_valid;
            end
        end
    end

    % Compte par groupe et condition (en excluant ceux 100% NaN)
    for ig = 1:numel(groups)
        gr = groups{ig};
        pids = fieldnames(contrib.(gr));
        if isempty(pids)
            for ic = 1:numel(conds)
                cond = conds{ic};
                if ~isfield(TABLES_AREA.N_Participant, gr), TABLES_AREA.N_Participant.(gr) = struct(); end
                if ~isfield(TABLES_AREA.N_Participant.(gr), pair), TABLES_AREA.N_Participant.(gr).(pair) = struct(); end
                TABLES_AREA.N_Participant.(gr).(pair).(cond) = 0;
            end
            continue;
        end

        M = false(numel(pids), numel(conds));
        for k = 1:numel(pids)
            M(k,:) = contrib.(gr).(pids{k});
        end

        keep = any(M, 2);
        M = M(keep, :);

        for ic = 1:numel(conds)
            cond = conds{ic};
            n_here = sum(M(:,ic));
            if ~isfield(TABLES_AREA.N_Participant, gr), TABLES_AREA.N_Participant.(gr) = struct(); end
            if ~isfield(TABLES_AREA.N_Participant.(gr), pair), TABLES_AREA.N_Participant.(gr).(pair) = struct(); end
            TABLES_AREA.N_Participant.(gr).(pair).(cond) = n_here;
        end
    end
end

%=========================== 6) Save master ===============================
mat_out = fullfile(stat_dir, 'ALL_TABLES_AREA_STRUCT.mat');
save(mat_out, 'TABLES_AREA', '-v7.3');
fprintf('DONE: %s\n', mat_out);

%% ============================== 2) PLOTS ================================
clc; close all;

% --- Chemins ---
root_coh = 'C:\Users\defsil00\Documents\Script\Results\Coherence';
all_dir  = fullfile(root_coh, 'ALL');                 % contains Coherence_*.mat
stat_dir = fullfile(root_coh, 'STATISTIQUE_AREA');    % new output root for AREA
csv_dir  = fullfile(stat_dir, 'CSV_CoherenceArea');
addpath(genpath('C:\Users\defsil00\Documents\Script'));
addpath(genpath('C:\Users\defsil00\Documents\Script'));

% --- Paramètres globaux ---
Conditions = {'Plat','Medium','High'};
Bandes     = {'Alpha','Beta','Gamma'};

% Libellés de bandes
band_ranges = struct();
band_ranges.Alpha = '8-12 Hz';
band_ranges.Beta  = '13-30 Hz';
band_ranges.Gamma = '31-60 Hz';

% Esthétique
col_left   = [0.00 0.60 0.00];  % vert
col_right  = [0.50 0.00 0.50];  % violet
alpha_pts  = 0.40;              % transparence des points
gray_levels = [0.85 0.85 0.85; 0.70 0.70 0.70; 0.55 0.55 0.55]; % Plat/Medium/High

% Dossier de sortie
fig_dir = fullfile(stat_dir, 'Fig_CoherenceArea_GroupedBars');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% Charger TABLES_AREA si nécessaire
if ~exist('TABLES_AREA','var')
    S = load(fullfile(stat_dir,'ALL_TABLES_AREA_STRUCT.mat'));
    TABLES_AREA = S.TABLES_AREA;
end

% Groupes d'âge (pour filtrage)
ParticipantGroup;

% Choisir le groupe à afficher
groupe_age = 'Adolescents'; % 'JeunesEnfants' | 'Enfants' | 'Adolescents' | 'Adultes'

% Récupérer les paires présentes (hors branches de synthèse)
all_pairs = setdiff(fieldnames(TABLES_AREA), {'N_cycle','N_Participant'});

% Lister toutes les sous-phases effectivement présentes
all_subphases = {};
for ip = 1:numel(all_pairs)
    p = all_pairs{ip};
    bnames = fieldnames(TABLES_AREA.(p));
    for ib = 1:numel(bnames)
        b = bnames{ib};
        spnames = fieldnames(TABLES_AREA.(p).(b));
        all_subphases = [all_subphases, spnames']; %#ok<AGROW>
    end
end
all_subphases = unique(all_subphases);

fprintf('=== Grouped bar plots for COHERENCE AREA (groupe: %s) ===\n', groupe_age);

for ib = 1:numel(Bandes)
    band = Bandes{ib};

    for isp = 1:numel(all_subphases)
        sp = all_subphases{isp};

        % Liste des paires qui ont bien les champs nécessaires
        pair_list = {};
        for ip = 1:numel(all_pairs)
            p = all_pairs{ip};
            if isfield(TABLES_AREA.(p), band) && isfield(TABLES_AREA.(p).(band), sp) ...
               && all(isfield(TABLES_AREA.(p).(band).(sp), {'left','right','mean'}))
                pair_list{end+1} = p; %#ok<AGROW>
            end
        end
        if isempty(pair_list)
            fprintf('  [!] Aucune paire pour subphase %s / bande %s\n', sp, band);
            continue;
        end

        % Matrices moyennes/SD : (lignes = conditions, colonnes = paires)
        M  = nan(numel(Conditions), numel(pair_list));
        SD = nan(numel(Conditions), numel(pair_list));

        % Stocker les tables par côté (pour les scatters)
        Tm_cell = cell(numel(pair_list),1);
        Tl_cell = cell(numel(pair_list),1);
        Tr_cell = cell(numel(pair_list),1);

        for ip = 1:numel(pair_list)
            p = pair_list{ip};
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
                vals = Tm.(c);
                M(ic, ip)  = mean(vals, 'omitnan');
                SD(ic, ip) = std(vals,  'omitnan'); % écart-type
            end
        end

        % ================== PLOT SUR AXE NUMÉRIQUE ==================
        x = 1:numel(pair_list);
        f = figure('Position',[100 120 1600 620]); hold on;  % un peu plus haut

        % Barres groupées (sur x numériques)
        bh = bar(x, M', 'grouped');   % M' = [nPairs x 3]
        for ic = 1:numel(Conditions)
            bh(ic).FaceColor = gray_levels(ic,:);
            bh(ic).EdgeColor = [0 0 0];
        end
        drawnow;

        % Centres des barres pour erreurs/points
        xC = arrayfun(@(b)b.XEndPoints, bh, 'UniformOutput', false);

        % Barres d'erreur + points de moyenne
        for ic = 1:numel(Conditions)
            xc = xC{ic};
            for ip = 1:numel(pair_list)
                mu = M(ic, ip); sd = SD(ic, ip);
                if ~isnan(mu)
                    plot(xc(ip), mu, 'ko', 'MarkerFaceColor','k', 'MarkerSize',6);
                    errorbar(xc(ip), mu, sd, 'k', 'LineWidth',1.2, 'CapSize',8);
                end
            end
        end

        % Points individuels par côté (avec masque NaN + jitter)
        jitter = 0.06;
        for ic = 1:numel(Conditions)
            xc = xC{ic};
            for ip = 1:numel(pair_list)
                % LEFT
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
                % RIGHT
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

        % ---------- Axe X : noms de paires uniquement ----------
        set(gca, 'XTick', x, 'XTickLabel', strrep(pair_list,'_','-'), ...
                 'TickLabelInterpreter','none');
        xlabel('Paires musculaires', 'FontSize',12, 'FontWeight','bold');
        ylabel('Coherence area (a.u.)', 'FontSize',12, 'FontWeight','bold');

        % ---------- Légende ----------
        hLeft  = scatter(nan, nan, 36, 'o', 'MarkerFaceColor',col_left,  'MarkerEdgeColor',col_left);
        hRight = scatter(nan, nan, 36, 'o', 'MarkerFaceColor',col_right, 'MarkerEdgeColor',col_right);
        legend([bh(1), bh(2), bh(3), hLeft, hRight], ...
               {'Plat','Medium','High','left','right'}, 'Location','northwest');

       % ---------- Panneau latéral avec N et C (PAR CONDITION) ----------
% Colonnes: PAIRE | N_Plat N_Med N_High | C_Plat C_Med C_High
wPair = 14;
wNs   = 5;  % largeur nombre N
wCs   = 6;  % largeur nombre C

fmtHeader = sprintf('%%-%ds  %%%ds %%%ds %%%ds   %%%ds %%%ds %%%ds', ...
                    wPair, wNs, wNs, wNs, wCs, wCs, wCs);
fmtRow    = sprintf('%%-%ds  %%%dd %%%dd %%%dd   %%%dd %%%dd %%%dd', ...
                    wPair, wNs, wNs, wNs, wCs, wCs, wCs);

info_lines = cell(1, numel(pair_list)+1);
info_lines{1} = sprintf(fmtHeader, 'PAIRE', 'N_P', 'N_M', 'N_H', 'C_P', 'C_M', 'C_H');

for ipx = 1:numel(pair_list)
    pp = strrep(string(pair_list{ipx}), '_', '-');
    Tm_pair = Tm_cell{ipx};

    % N participants PAR CONDITION (sur la table 'mean' filtrée au groupe)
    Np = [0 0 0];  % [Plat Medium High]
    if ~isempty(Tm_pair)
        for ic = 1:numel(Conditions)
            c = Conditions{ic};
            if ismember(c, Tm_pair.Properties.VariableNames)
                Np(ic) = sum(~isnan(Tm_pair.(c)));
            end
        end
    end

    % C cycles PAR CONDITION (from TABLES_AREA.N_cycle)
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

% --- Réserver un peu plus d'espace pour le panneau à droite
ax = gca;
ax.Position = [0.08 0.18 0.56 0.72];   % <- largeur de l'axe réduite (0.62 -> 0.56)

% --- Titre (inchangé)
t = title(sprintf('%s (%s) - %s - %s', band, band_ranges.(band), sp, groupe_age), ...
          'FontSize',14,'FontWeight','bold');
t.Units = 'normalized'; t.Position(2) = 1.03;

% --- Annotation (panneau N / C) : élargie
annotation(f, 'textbox', [0.68 0.18 0.32 0.72], ...   % <- largeur 0.32 (au lieu de 0.26)
           'String', strjoin(info_lines, newline), ...
           'Interpreter','none', ...
           'FontName','Consolas', 'FontSize',10, ...
           'HorizontalAlignment','left', ...
           'VerticalAlignment','top', ...
           'EdgeColor',[0.85 0.85 0.85], 'BackgroundColor',[1 1 1], ...
           'FitBoxToText','off', ...
           'Margin', 6);                               % un peu d'air mais pas trop

        % ---------- Export propre ----------
        figname = sprintf('GroupedBars_AREA_%s_%s_%s.png', sp, band, groupe_age);
        exportgraphics(f, fullfile(fig_dir, figname), 'Resolution', 300);
        fprintf('  OK FIG: %s / %s -> %s\n', sp, band, figname);
        close(f);
    end
end

fprintf('All done.\n');

%% =========================== LOCAL FUNCTIONS ============================
function grp = resolveAge(pid, Group)
% Return age group name for a given participant id (string)
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
% Retourne L_target_eff_<pair> pour CETTE condition si ≥1 cycle retenu (G ou D), sinon 0.

    n = 0;
    if ~isstruct(DATA) || ~isfield(DATA, cond), return; end

    parts = strsplit(pair,'_');
    if numel(parts) ~= 2, return; end
    m1 = parts{1}; m2 = parts{2};

    % 1) Y a-t-il au moins un cycle retenu sur ce cond pour cette paire ?
    has_cycles = local_side_has_cycles(DATA.(cond), 'left',  m1, m2) ...
              || local_side_has_cycles(DATA.(cond), 'right', m1, m2);
    if ~has_cycles
        n = 0;   % condition vide -> quota = 0 pour ce participant
        return;
    end

    % 2) Lire le quota d’étude dans la méta de CETTE condition
    if isfield(DATA.(cond), 'meta') && isstruct(DATA.(cond).meta)
        meta = DATA.(cond).meta;
        f1 = ['L_target_eff_' pair];
        f2 = ['L_target_'     pair];   % fallback si anciens fichiers
        if isfield(meta, f1) && isnumeric(meta.(f1)) && isfinite(meta.(f1))
            n = double(meta.(f1)); return;
        elseif isfield(meta, f2) && isnumeric(meta.(f2)) && isfinite(meta.(f2))
            n = double(meta.(f2)); return;
        end
    end

    % Si pas trouvable: applique ta règle "non trouvé => 0"
    n = 0;
end

function tf = local_side_has_cycles(node_cond, side, m1, m2)
% True si ce côté a ≥1 cycle retenu (on lit Ncycle_<m1>_<m2> ou à défaut L_<m1>_<m2>)
    tf = false;
    if ~isfield(node_cond, side), return; end
    S = node_cond.(side);

    fnN = ['Ncycle_' m1 '_' m2];
    if isfield(S, fnN) && isnumeric(S.(fnN)) && isfinite(S.(fnN)) && S.(fnN) > 0
        tf = true; return;
    end
    fnL = ['L_' m1 '_' m2];
    if isfield(S, fnL) && isnumeric(S.(fnL)) && isfinite(S.(fnL)) && S.(fnL) > 0
        tf = true; return;
    end
end