%% VISUALISATION DU NIVEAU D'ACTIVATION MUSCULAIRE PAR CONDITION ET SOUS-PHASE
%
% OBJECTIF :
%   Ce script calcule et visualise le profil d'activation musculaire moyen
%   sur le cycle de marche pour un groupe d'âge donné. Il produit, pour chaque
%   muscle, une courbe moyenne normalisée par condition (Plat, Medium, High)
%   avec ses bandes d'écart-type, superposée aux événements de marche clés
%   (Toe-Off principal, Toe-Off et Heel Strike controlatéraux).
%
%   PIPELINE :
%     1. Chargement des cycles moyens pré-calculés depuis *_MATRIX.mat
%        (structure CYCLES_MOYENS, jambes gauche et droite confondues)
%     2. Calcul des cycles moyens par participant → moyenne de groupe
%     3. Normalisation par le maximum de la condition Plat (par muscle)
%     4. Calcul des événements de marche moyens du groupe :
%        - Toe-Off principal (depuis CYCLES_TOEOFF dans *_MATRIX.mat)
%        - Toe-Off et Heel Strike controlatéraux (depuis les fichiers
%          Coherence_*.mat, champs GaitEvents_*)
%     5. Génération d'une figure à 8 sous-graphes (un par muscle) avec
%        courbes par condition, bandes d'écart-type et marqueurs d'événements
%
% ENTRÉES :
%   - *_MATRIX.mat   : Fichiers pré-calculés par participant (dossier FILTERED)
%                      contenant :
%       * CYCLES_MOYENS.(pid).(cond).(side).(muscle) : Cycle moyen normalisé
%       * CYCLES_TOEOFF.(pid).(cond).(side).mean_percentage : % du TO
%   - Coherence_<PID>.mat : Fichiers du script 7, pour les événements
%                            controlatéraux (GaitEvents_*)
%   - ParticipantGroup.m  : Définit la struct Group (listes de participants)
%
% SORTIES :
%   - Cycle_Normalise_AllConditions_AllMuscles.png :
%     Figure à 8 sous-graphes (un par muscle), avec pour chaque condition :
%       * Courbe de l'activation normalisée moyenne ± écart-type (bande)
%       * Ligne verticale pointillée au Toe-Off principal (--) 
%       * Ligne verticale pointillée au Toe-Off controlatéral (:)
%       * Ligne verticale pointillée au Heel Strike controlatéral (-.)
%     Sauvegardée dans : Results/Coherence/Visualisation_Coherence/<groupe>/
%
% NORMALISATION :
%   Chaque muscle est normalisé par le maximum de son cycle moyen de groupe
%   en condition Plat. Ceci permet de comparer l'amplitude relative entre
%   conditions pour un même muscle (1.0 = pic d'activation sur Plat).
%
% NOTES :
%   - Les deux jambes (gauche et droite) sont moyennées ensemble par participant
%     avant d'être incluses dans la moyenne de groupe.
%   - Les événements sont légèrement décalés horizontalement (jitter) pour
%     désuperposer visuellement les marqueurs des 3 conditions.
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CONFIGURATION ====================================

% Répertoire contenant les fichiers *_MATRIX.mat
cd('C:\Users\defsil00\Documents\Script\FILTERED');
addpath(genpath('C:\Users\defsil00\Documents\Script'));

% Groupe d'âge à analyser
groupe_a_etudier = 'Enfants'; % 'JeunesEnfants' | 'Enfants' | 'Adolescents' | 'Adultes'

% Chargement de la définition des groupes
% → Crée la variable struct Group avec les listes de participants
ParticipantGroup;
Participant = Group.(groupe_a_etudier); % Liste des IDs du groupe sélectionné

fprintf('\n=== Groupe sélectionné : %s ===\n', groupe_a_etudier);
fprintf('Participants inclus : %s\n\n', strjoin(Participant, ', '));

%% ===================== PARAMÈTRES =======================================

Conditions = {'Plat','Medium','High'};
Muscles    = {'EMG_TAprox','EMG_TAdist','EMG_SOL','EMG_GM', ...
              'EMG_VL','EMG_RF','EMG_ST','EMG_GMED'};
Jambes     = {'left','right'};

%% ===================== INITIALISATIONS ==================================

% mean_cycles_combined.(cond).(muscle) : matrice [nParticipants × 101 pts]
% Chaque ligne = cycle moyen d'un participant (jambes confondues)
mean_cycles_combined = struct();

% participant_tracking.(cond).(muscle) : cell array des PIDs correspondant
% aux lignes de mean_cycles_combined (utile pour le débogage)
participant_tracking = struct();

%% ===================== CHARGEMENT DES CYCLES MOYENS =====================

for iP = 1:length(Participant)
    filename = [Participant{iP} '_MATRIX.mat'];
    if ~exist(filename,'file')
        fprintf('Fichier manquant: %s\n', filename);
        continue;
    end
    metadonnees = load(filename);

    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        for iM = 1:length(Muscles)
            muscle = Muscles{iM};

            % Accumulation des cycles valides des deux jambes
            valid_cycles = [];
            for j = 1:length(Jambes)
                jambe = Jambes{j};

                % Vérification de l'existence du champ avant accès
                if isfield(metadonnees,'CYCLES_MOYENS') && ...
                   isfield(metadonnees.CYCLES_MOYENS, Participant{iP}) && ...
                   isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}), condition) && ...
                   isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition), jambe) && ...
                   isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe), muscle)

                    cycle_mean = metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe).(muscle);

                    if all(isfinite(cycle_mean))
                        valid_cycles = [valid_cycles; cycle_mean]; 
                        fprintf('Cycle valide: %s - %s - %s - %s\n', Participant{iP}, condition, jambe, muscle);
                    else
                        fprintf('Cycle NaN ignoré: %s - %s - %s - %s\n', Participant{iP}, condition, jambe, muscle);
                    end
                else
                    fprintf('Muscle absent: %s - %s - %s - %s\n', muscle, Participant{iP}, condition, jambe);
                end
            end

            % Calcul du cycle moyen du participant (jambes confondues)
            if ~isempty(valid_cycles)
                participant_mean_cycle = mean(valid_cycles, 1);
                if all(isfinite(participant_mean_cycle))

                    % Initialisation des champs si nécessaire
                    if ~isfield(mean_cycles_combined, condition)
                        mean_cycles_combined.(condition) = struct();
                        participant_tracking.(condition) = struct();
                    end
                    if ~isfield(mean_cycles_combined.(condition), muscle)
                        mean_cycles_combined.(condition).(muscle) = [];
                        participant_tracking.(condition).(muscle) = {};
                    end

                    % Empilement du cycle moyen du participant
                    mean_cycles_combined.(condition).(muscle) = ...
                        [mean_cycles_combined.(condition).(muscle); participant_mean_cycle];
                    participant_tracking.(condition).(muscle){end+1} = Participant{iP};
                    fprintf('Ajouté : %s -> %s - %s\n', Participant{iP}, condition, muscle);
                else
                    fprintf('Cycle moyen NaN rejeté: %s - %s - %s\n', Participant{iP}, condition, muscle);
                end
            else
                fprintf('Aucun cycle valide: %s - %s - %s\n', Participant{iP}, condition, muscle);
            end
        end
    end
end

% ===================== MOYENNES DE GROUPE ================================
% Calcul de la moyenne sur tous les participants pour chaque condition/muscle

global_mean_cycles = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    global_mean_cycles.(condition) = struct();

    for iM = 1:length(Muscles)
        muscle = Muscles{iM};
        if isfield(mean_cycles_combined, condition) && ...
           isfield(mean_cycles_combined.(condition), muscle) && ...
           ~isempty(mean_cycles_combined.(condition).(muscle))

            tmp = mean(mean_cycles_combined.(condition).(muscle), 1);
            if ~any(isnan(tmp))
                global_mean_cycles.(condition).(muscle) = tmp;
                fprintf('Moyenne globale calculée: %s - %s\n', condition, muscle);
            end
        end
    end
end

% ===================== NORMALISATION PAR LE MAX DE PLAT =================
% Normalise chaque muscle par le maximum de son cycle moyen de groupe en
% condition Plat → l'amplitude de Plat vaut 1.0 pour chaque muscle

normalized_mean_cycles = struct();
for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    if isfield(global_mean_cycles,'Plat') && isfield(global_mean_cycles.Plat, muscle)
        max_flat = max(global_mean_cycles.Plat.(muscle));

        if ~isnan(max_flat) && max_flat > 0
            for iC = 1:length(Conditions)
                condition = Conditions{iC};
                if isfield(global_mean_cycles, condition) && ...
                   isfield(global_mean_cycles.(condition), muscle)
                    if ~isfield(normalized_mean_cycles, condition)
                        normalized_mean_cycles.(condition) = struct();
                    end
                    normalized_mean_cycles.(condition).(muscle) = ...
                        global_mean_cycles.(condition).(muscle) / max_flat;
                    fprintf('Normalisé: %s - %s\n', condition, muscle);
                end
            end
        else
            fprintf('Max Plat invalide pour %s (%.3f) → pas de normalisation\n', muscle, max_flat);
        end
    else
        fprintf('Condition Plat manquante pour %s → pas de normalisation\n', muscle);
    end
end

% ===================== NORMALISATION PAR PARTICIPANT (pour figures indiv.) =
% Version de la normalisation maintenue à l'échelle du participant
% (utile pour des figures individuelles, non utilisée dans la figure de groupe)

normalized_participant_cycles = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    if isfield(mean_cycles_combined, condition)
        normalized_participant_cycles.(condition) = struct();
        for iM = 1:length(Muscles)
            muscle = Muscles{iM};
            if isfield(mean_cycles_combined.(condition), muscle) && ...
               isfield(global_mean_cycles,'Plat') && isfield(global_mean_cycles.Plat, muscle)
                max_flat = max(global_mean_cycles.Plat.(muscle));
                if ~isnan(max_flat) && max_flat > 0
                    normalized_participant_cycles.(condition).(muscle) = ...
                        mean_cycles_combined.(condition).(muscle) / max_flat;
                end
            end
        end
    end
end

% ===================== CALCUL DU TOE-OFF MOYEN DU GROUPE ================
% Pour chaque condition, calcule le pourcentage moyen du cycle auquel
% survient le Toe-Off principal (depuis CYCLES_TOEOFF dans *_MATRIX.mat).
% Les jambes gauche et droite sont moyennées ; une seule jambe disponible
% est acceptée (robuste aux données unilatérales).

toeoff_moyen = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    toeoff_condition = [];

    for iP = 1:length(Participant)
        participant = Participant{iP};
        filename = [participant '_MATRIX.mat'];
        if ~exist(filename,'file'), continue; end
        data = load(filename);

        left_toeoff  = NaN;
        right_toeoff = NaN;

        try
            if isfield(data,'CYCLES_TOEOFF') && ...
               isfield(data.CYCLES_TOEOFF, participant) && ...
               isfield(data.CYCLES_TOEOFF.(participant), condition)

                % Lecture du TO gauche si disponible
                if isfield(data.CYCLES_TOEOFF.(participant).(condition),'left') && ...
                   isfield(data.CYCLES_TOEOFF.(participant).(condition).left,'mean_percentage')
                    left_toeoff = data.CYCLES_TOEOFF.(participant).(condition).left.mean_percentage;
                end
                % Lecture du TO droit si disponible
                if isfield(data.CYCLES_TOEOFF.(participant).(condition),'right') && ...
                   isfield(data.CYCLES_TOEOFF.(participant).(condition).right,'mean_percentage')
                    right_toeoff = data.CYCLES_TOEOFF.(participant).(condition).right.mean_percentage;
                end
            end
        catch
            % Participant sans données de TO : ignoré silencieusement
        end

        % Accepte les valeurs unilatérales (nanmean ignore les NaN)
        if any(isfinite([left_toeoff right_toeoff]))
            toeoff_condition(end+1) = mean([left_toeoff right_toeoff], 'omitnan');
        end
    end

    if ~isempty(toeoff_condition)
        toeoff_moyen.(condition) = mean(toeoff_condition, 'omitnan');
        fprintf('Toe-Off moyen %s = %.2f%% (N=%d)\n', condition, toeoff_moyen.(condition), numel(toeoff_condition));
    else
        toeoff_moyen.(condition) = NaN;
        fprintf('Toe-Off moyen %s = NaN (aucun participant valide)\n', condition);
    end
end

% ===================== CALCUL DES ÉVÉNEMENTS CONTROLATÉRAUX =============
% Récupère les TO et HS controlatéraux moyens depuis les champs GaitEvents_*
% des fichiers Coherence_<PID>.mat (calculés dans le script 7).

coh_all_dir           = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
oppo_toeoff_moyen     = struct();
oppo_heelstrike_moyen = struct();

for iC = 1:length(Conditions)
    condition    = Conditions{iC};
    vals_oppo_to = [];
    vals_oppo_hs = [];

    for iP = 1:length(Participant)
        pid   = Participant{iP};
        f_coh = fullfile(coh_all_dir, ['Coherence_' pid '.mat']);
        if ~exist(f_coh,'file'), continue; end
        S = load(f_coh);
        if ~isfield(S,'DATA') || ~isfield(S.DATA, condition), continue; end

        sides        = {'left','right'};
        oppo_to_lr   = nan(1,2);
        oppo_hs_lr   = nan(1,2);

        for s = 1:2
            sd = sides{s};
            if ~isfield(S.DATA.(condition), sd), continue; end

            % Recherche du premier champ GaitEvents_* disponible
            fns = fieldnames(S.DATA.(condition).(sd));
            idx = find(startsWith(fns,'GaitEvents_'), 1, 'first');
            if isempty(idx)
                % Fallback vers le champ de référence TAprox_TAdist
                cand = 'GaitEvents_TAprox_TAdist';
                if isfield(S.DATA.(condition).(sd), cand)
                    ge = S.DATA.(condition).(sd).(cand);
                else
                    continue;
                end
            else
                ge = S.DATA.(condition).(sd).(fns{idx});
            end

            % Extraction des événements controlatéraux
            if isstruct(ge)
                if isfield(ge,'opposite_toeoff') && isfinite(ge.opposite_toeoff)
                    oppo_to_lr(s) = ge.opposite_toeoff;
                end
                if isfield(ge,'opposite_heelstrike') && isfinite(ge.opposite_heelstrike)
                    oppo_hs_lr(s) = ge.opposite_heelstrike;
                end
            end
        end

        % Moyenne bilatérale (accepte les valeurs unilatérales)
        if any(isfinite(oppo_to_lr))
            vals_oppo_to(end+1) = mean(oppo_to_lr,'omitnan'); 
        end
        if any(isfinite(oppo_hs_lr))
            vals_oppo_hs(end+1) = mean(oppo_hs_lr,'omitnan'); 
        end
    end

    % Calcul des moyennes de groupe
    if ~isempty(vals_oppo_to)
        oppo_toeoff_moyen.(condition) = mean(vals_oppo_to,'omitnan');
    else
        oppo_toeoff_moyen.(condition) = NaN;
    end
    if ~isempty(vals_oppo_hs)
        oppo_heelstrike_moyen.(condition) = mean(vals_oppo_hs,'omitnan');
    else
        oppo_heelstrike_moyen.(condition) = NaN;
    end

    fprintf('[%s] oppo_TO=%.2f%% | oppo_HS=%.2f%%\n', condition, ...
        oppo_toeoff_moyen.(condition), oppo_heelstrike_moyen.(condition));
end

% ===================== FIGURE COMBINÉE (3 CONDITIONS × 8 MUSCLES) =======

% Dossier de sortie propre au groupe
base_output_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\Visualisation_Coherence';
output_dir = fullfile(base_output_dir, groupe_a_etudier);
if ~exist(output_dir,'dir'), mkdir(output_dir); end

% Couleurs pour les 3 conditions
colors = {'b','g','r'}; % Plat = bleu, Medium = vert, High = rouge

figure('Name','Cycle - 3 Conditions','Position',[100 100 1800 900]);

for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    subplot(2,4,iM); hold on;
    title(strrep(muscle,'EMG_',''), 'Interpreter','none'); % Supprime le préfixe 'EMG_'
    xlabel('Cycle (%)'); ylabel('Activation normalisée');
    ylim([0 2]);

    h_legend = [];

    % --- Tracé des courbes d'activation par condition ---
    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        if isfield(normalized_mean_cycles, condition) && ...
           isfield(normalized_mean_cycles.(condition), muscle)

            cycle = normalized_mean_cycles.(condition).(muscle);

            % Calcul de l'écart-type normalisé pour la bande d'incertitude
            if isfield(global_mean_cycles,'Plat') && isfield(global_mean_cycles.Plat, muscle)
                max_flat = max(global_mean_cycles.Plat.(muscle));
                if ~isnan(max_flat) && max_flat > 0
                    vals = mean_cycles_combined.(condition).(muscle);
                    if isempty(vals)
                        std_cycle = zeros(size(cycle));
                    else
                        std_cycle = std(vals, 0, 1, 'omitnan') / max_flat;
                    end
                else
                    std_cycle = zeros(size(cycle));
                end
            else
                std_cycle = zeros(size(cycle));
            end

            % Courbe moyenne
            h = plot(cycle, 'LineWidth',2, 'Color',colors{iC});

            % Bande ± écart-type (transparente)
            fill([1:numel(cycle), fliplr(1:numel(cycle))], ...
                 [cycle + std_cycle, fliplr(cycle - std_cycle)], ...
                 colors{iC}, 'FaceAlpha',0.2, 'EdgeColor','none');
            h_legend(end+1) = h; 
        end
    end

    % --- Tracé des événements de marche ---
    % Détermination de la longueur de référence du vecteur cycle (101 points)
    npts_ref = [];
    for c = 1:length(Conditions)
        cond = Conditions{c};
        if isfield(normalized_mean_cycles,cond) && isfield(normalized_mean_cycles.(cond), muscle)
            npts_ref = numel(normalized_mean_cycles.(cond).(muscle));
            break
        end
    end
    if isempty(npts_ref), npts_ref = 101; end
    yl = ylim;

    % Jitter horizontal pour désuperposer visuellement les traits des 3 conditions
    event_jitter_pct = [-0.6, 0, +0.6]; % [Plat, Medium, High] en % du cycle

    % Conversion % du cycle → index sur le vecteur (avec jitter)
    pct2idx = @(pct,j) max(1, min(npts_ref, round((pct + j)/100 * npts_ref)));

    for c = 1:length(Conditions)
        cond = Conditions{c};
        col  = colors{c};
        j    = event_jitter_pct(c);

        % Toe-Off principal (ligne tiretée --) + marqueur triangle bas
        if isfield(toeoff_moyen,cond) && ~isnan(toeoff_moyen.(cond))
            x = pct2idx(toeoff_moyen.(cond), j);
            plot([x x], yl, '--', 'Color',col, 'LineWidth',1.5);
            scatter(x, yl(2), 30, col, 'v', 'filled', 'MarkerEdgeColor','w');
            text(x, yl(2), sprintf(' %.1f%%', toeoff_moyen.(cond)), ...
                'Color',col, 'FontSize',8, 'VerticalAlignment','top', 'FontWeight','bold');
        end

        % Toe-Off controlatéral (ligne pointillée :) + marqueur triangle haut
        if isfield(oppo_toeoff_moyen,cond) && ~isnan(oppo_toeoff_moyen.(cond))
            x = pct2idx(oppo_toeoff_moyen.(cond), j);
            plot([x x], yl, ':', 'Color',col, 'LineWidth',1.5);
            scatter(x, yl(2), 30, col, '^', 'filled', 'MarkerEdgeColor','w');
            text(x, yl(2), sprintf(' %.1f%%', oppo_toeoff_moyen.(cond)), ...
                'Color',col, 'FontSize',8, 'VerticalAlignment','bottom', 'FontWeight','bold');
        end

        % Heel Strike controlatéral (ligne dash-point -.) + marqueur carré
        if isfield(oppo_heelstrike_moyen,cond) && ~isnan(oppo_heelstrike_moyen.(cond))
            x = pct2idx(oppo_heelstrike_moyen.(cond), j);
            plot([x x], yl, '-.', 'Color',col, 'LineWidth',1.5);
            scatter(x, yl(2), 30, col, 's', 'filled', 'MarkerEdgeColor','w');
            text(x, yl(2), sprintf(' %.1f%%', oppo_heelstrike_moyen.(cond)), ...
                'Color',col, 'FontSize',8, 'VerticalAlignment','middle', 'FontWeight','bold');
        end
    end

    % --- Formatage de l'axe X (graduation 0 à 100% par pas de 10%) ---
    xt = round(linspace(1, npts_ref, 11));
    xticks(xt);
    xticklabels(compose('%d%%', 0:10:100));
    xlim([1 npts_ref]);

    grid on;

    % Légende uniquement sur le premier sous-graphe
    if iM == 1 && ~isempty(h_legend)
        legend(h_legend, Conditions, 'Location','best');
    end
end

%% ===================== SAUVEGARDE =======================================
saveas(gcf, fullfile(output_dir, 'Cycle_Normalise_AllConditions_AllMuscles.png'));
fprintf('Figure sauvegardée : Cycle_Normalise_AllConditions_AllMuscles.png\n');

fprintf('\nTerminé.\n');