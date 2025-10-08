%% VISUALISATION ACTIVATION MUSCULAIRE (corrige)
clc; clear; close all;

% Paths
cd('C:\Users\defsil00\Documents\Script\FILTERED');
addpath(genpath('C:\Users\defsil00\Documents\Script'));

% === CONFIGURATION DU GROUPE A ETUDIER ===
groupe_a_etudier = 'Enfants';  % 'JeunesEnfants' | 'Enfants' | 'Adolescents' | 'Adultes'

% Charger la structure Group definie dans ParticipantGroup.m
ParticipantGroup;

% Liste des participants du groupe
Participant = Group.(groupe_a_etudier);

% Affichage
fprintf('\n=== Groupe selectionne : %s ===\n', groupe_a_etudier);
fprintf('Participants inclus : %s\n\n', strjoin(Participant, ', '));

Conditions = {'Plat','Medium','High'};
Muscles    = {'EMG_TAprox','EMG_TAdist','EMG_SOL','EMG_GM','EMG_VL','EMG_RF','EMG_ST','EMG_GMED'};
Jambes     = {'left','right'};

% Stockage
mean_cycles_combined  = struct();   % cycles moyens par participant, puis empiles
participant_tracking  = struct();   % noms des participants alignes aux lignes

% ====== Charger cycles moyens par participant / condition / muscle ======
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

            valid_cycles = [];
            for j = 1:length(Jambes)
                jambe = Jambes{j};
                if isfield(metadonnees,'CYCLES_MOYENS') && ...
                   isfield(metadonnees.CYCLES_MOYENS, Participant{iP}) && ...
                   isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}), condition) && ...
                   isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition), jambe) && ...
                   isfield(metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe), muscle)

                    cycle_mean = metadonnees.CYCLES_MOYENS.(Participant{iP}).(condition).(jambe).(muscle);
                    if all(isfinite(cycle_mean))
                        valid_cycles = [valid_cycles; cycle_mean]; %#ok<AGROW>
                        fprintf('Cycle valide: %s - %s - %s - %s\n', Participant{iP}, condition, jambe, muscle);
                    else
                        fprintf('Cycle NaN ignore: %s - %s - %s - %s\n', Participant{iP}, condition, jambe, muscle);
                    end
                else
                    fprintf('Muscle absent: %s - %s - %s - %s\n', muscle, Participant{iP}, condition, jambe);
                end
            end

            if ~isempty(valid_cycles)
                participant_mean_cycle = mean(valid_cycles,1);
                if all(isfinite(participant_mean_cycle))
                    if ~isfield(mean_cycles_combined, condition)
                        mean_cycles_combined.(condition) = struct();
                        participant_tracking.(condition) = struct();
                    end
                    if ~isfield(mean_cycles_combined.(condition), muscle)
                        mean_cycles_combined.(condition).(muscle) = [];
                        participant_tracking.(condition).(muscle) = {};
                    end

                    mean_cycles_combined.(condition).(muscle) = ...
                        [mean_cycles_combined.(condition).(muscle); participant_mean_cycle]; %#ok<AGROW>
                    participant_tracking.(condition).(muscle){end+1} = Participant{iP}; %#ok<AGROW>

                    fprintf('Ajoute participant %s -> %s - %s\n', Participant{iP}, condition, muscle);
                else
                    fprintf('Cycle moyen NaN rejete: %s - %s - %s\n', Participant{iP}, condition, muscle);
                end
            else
                fprintf('Aucun cycle valide: %s - %s - %s\n', Participant{iP}, condition, muscle);
            end
        end
    end
end

% ====== Moyennes globales (toutes jambes confondues) ======
global_mean_cycles = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    global_mean_cycles.(condition) = struct();

    for iM = 1:length(Muscles)
        muscle = Muscles{iM};
        if isfield(mean_cycles_combined, condition) && isfield(mean_cycles_combined.(condition), muscle)
            if ~isempty(mean_cycles_combined.(condition).(muscle))
                tmp = mean(mean_cycles_combined.(condition).(muscle),1);
                if ~any(isnan(tmp))
                    global_mean_cycles.(condition).(muscle) = tmp;
                    fprintf('Global mean: %s - %s\n', condition, muscle);
                end
            end
        end
    end
end

% ====== Normalisation par le max de Plat (par muscle) ======
normalized_mean_cycles = struct();
for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    if isfield(global_mean_cycles,'Plat') && isfield(global_mean_cycles.Plat, muscle)
        max_flat = max(global_mean_cycles.Plat.(muscle));
        if ~isnan(max_flat) && max_flat > 0
            for iC = 1:length(Conditions)
                condition = Conditions{iC};
                if isfield(global_mean_cycles, condition) && isfield(global_mean_cycles.(condition), muscle)
                    if ~isfield(normalized_mean_cycles, condition)
                        normalized_mean_cycles.(condition) = struct();
                    end
                    normalized_mean_cycles.(condition).(muscle) = global_mean_cycles.(condition).(muscle) / max_flat;
                    fprintf('Normalise: %s - %s\n', condition, muscle);
                end
            end
        else
            fprintf('Max Plat invalide pour %s (%.3f)\n', muscle, max_flat);
        end
    else
        fprintf('Plat manquant pour %s -> pas de normalisation\n', muscle);
    end
end

% ====== Normalisation par participant (optionnel pour figures indiv.) ======
normalized_participant_cycles = struct();
for iC = 1:length(Conditions)
    condition = Conditions{iC};
    if isfield(mean_cycles_combined, condition)
        normalized_participant_cycles.(condition) = struct();
        for iM = 1:length(Muscles)
            muscle = Muscles{iM};
            if isfield(mean_cycles_combined.(condition), muscle)
                if isfield(global_mean_cycles,'Plat') && isfield(global_mean_cycles.Plat, muscle)
                    max_flat = max(global_mean_cycles.Plat.(muscle));
                    if ~isnan(max_flat) && max_flat > 0
                        normalized_participant_cycles.(condition).(muscle) = ...
                            mean_cycles_combined.(condition).(muscle) / max_flat;
                    end
                end
            end
        end
    end
end

% ====== TOE-OFF groupe robuste (unilateral accepte) ======
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
            if isfield(data,'CYCLES_TOEOFF') && isfield(data.CYCLES_TOEOFF, participant) && ...
               isfield(data.CYCLES_TOEOFF.(participant), condition)

                if isfield(data.CYCLES_TOEOFF.(participant).(condition),'left') && ...
                   isfield(data.CYCLES_TOEOFF.(participant).(condition).left,'mean_percentage')
                    left_toeoff = data.CYCLES_TOEOFF.(participant).(condition).left.mean_percentage;
                end
                if isfield(data.CYCLES_TOEOFF.(participant).(condition),'right') && ...
                   isfield(data.CYCLES_TOEOFF.(participant).(condition).right,'mean_percentage')
                    right_toeoff = data.CYCLES_TOEOFF.(participant).(condition).right.mean_percentage;
                end
            end
        catch
            % ignore
        end

        if any(isfinite([left_toeoff right_toeoff]))
            toeoff_condition(end+1) = mean([left_toeoff right_toeoff], 'omitnan'); %#ok<AGROW>
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

% ====== Opposite events groupe (unilateral accepte) ======
coh_all_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\ALL';
oppo_toeoff_moyen     = struct();
oppo_heelstrike_moyen = struct();

for iC = 1:length(Conditions)
    condition = Conditions{iC};
    vals_oppo_to = [];
    vals_oppo_hs = [];

    for iP = 1:length(Participant)
        pid  = Participant{iP};
        f_coh = fullfile(coh_all_dir, ['Coherence_' pid '.mat']);
        if ~exist(f_coh,'file'), continue; end
        S = load(f_coh);
        if ~isfield(S,'DATA') || ~isfield(S.DATA,condition), continue; end

        sides = {'left','right'};
        oppo_to_lr = nan(1,2);
        oppo_hs_lr = nan(1,2);

        for s = 1:2
            sd = sides{s};
            if ~isfield(S.DATA.(condition), sd), continue; end

            % trouver un champ GaitEvents_*
            fns = fieldnames(S.DATA.(condition).(sd));
            idx = find(startsWith(fns,'GaitEvents_'),1,'first');
            if isempty(idx)
                cand = 'GaitEvents_TAprox_TAdist';
                if isfield(S.DATA.(condition).(sd), cand)
                    ge = S.DATA.(condition).(sd).(cand);
                else
                    continue;
                end
            else
                ge = S.DATA.(condition).(sd).(fns{idx});
            end

            if isstruct(ge)
                if isfield(ge,'opposite_toeoff') && isfinite(ge.opposite_toeoff)
                    oppo_to_lr(s) = ge.opposite_toeoff;
                end
                if isfield(ge,'opposite_heelstrike') && isfinite(ge.opposite_heelstrike)
                    oppo_hs_lr(s) = ge.opposite_heelstrike;
                end
            end
        end

        if any(isfinite(oppo_to_lr))
            vals_oppo_to(end+1) = mean(oppo_to_lr,'omitnan'); %#ok<AGROW>
        end
        if any(isfinite(oppo_hs_lr))
            vals_oppo_hs(end+1) = mean(oppo_hs_lr,'omitnan'); %#ok<AGROW>
        end
    end

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

% ================= FIGURES =================

% Output dir par groupe
base_output_dir = 'C:\Users\defsil00\Documents\Script\Results\Coherence\Visualisation_Coherence';
output_dir = fullfile(base_output_dir, groupe_a_etudier);
if ~exist(output_dir,'dir'), mkdir(output_dir); end

% Couleurs pour conditions (courbes globales)
colors = {'b','g','r'}; % Plat, Medium, High

% ====== FIGURE COMBINÉE AVEC LES 3 CONDITIONS ======
figure('Name','Cycle - 3 Conditions','Position',[100 100 1800 900]);
for iM = 1:length(Muscles)
    muscle = Muscles{iM};
    subplot(2,4,iM); hold on;
    title(strrep(muscle,'EMG_',''), 'Interpreter','none');
    xlabel('Cycle (%)'); ylabel('Activation normalisee');
    ylim([0 2]);

    h_legend = [];
    % tracer les courbes par condition si dispo
    for iC = 1:length(Conditions)
        condition = Conditions{iC};
        if isfield(normalized_mean_cycles,condition) && isfield(normalized_mean_cycles.(condition), muscle)
            cycle = normalized_mean_cycles.(condition).(muscle);
            if isfield(global_mean_cycles,'Plat') && isfield(global_mean_cycles.Plat, muscle)
                max_flat = max(global_mean_cycles.Plat.(muscle));
                if ~isnan(max_flat) && max_flat>0
                    % ==== REMPLACE nanstd PAR std(...,'omitnan') ====
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
            h = plot(cycle,'LineWidth',2,'Color',colors{iC});
            fill([1:numel(cycle), fliplr(1:numel(cycle))], ...
                 [cycle + std_cycle, fliplr(cycle - std_cycle)], ...
                 colors{iC}, 'FaceAlpha',0.2,'EdgeColor','none');
            h_legend(end+1) = h; %#ok<AGROW>
        end
    end

    % === Tracer les evenements POUR TOUTES LES CONDITIONS (hors if) ===
    % longueur de reference pour convertir % -> index
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

    % Décalage horizontal léger pour désuperposer visuellement les traits
    event_jitter_pct = [-0.6, 0, +0.6];  % Plat, Medium, High

    % Helper: % -> index avec jitter
    pct2idx = @(pct,j) max(1, min(npts_ref, round((pct + j)/100 * npts_ref)));

    for c = 1:length(Conditions)
        cond = Conditions{c};
        col  = colors{c};
        j    = event_jitter_pct(c);

        % Toe-off
        if isfield(toeoff_moyen,cond) && ~isnan(toeoff_moyen.(cond))
            x = pct2idx(toeoff_moyen.(cond), j);
            plot([x x], yl, '--', 'Color', col, 'LineWidth', 1.5);
            scatter(x, yl(2), 30, col, 'v', 'filled', 'MarkerEdgeColor','w'); % marqueur (optionnel)
            text(x, yl(2), sprintf(' %.1f%%', toeoff_moyen.(cond)), ...
                'Color', col, 'FontSize', 8, 'VerticalAlignment','top', 'FontWeight','bold');
        end

        % Opposite toe-off
        if isfield(oppo_toeoff_moyen,cond) && ~isnan(oppo_toeoff_moyen.(cond))
            x = pct2idx(oppo_toeoff_moyen.(cond), j);
            plot([x x], yl, ':', 'Color', col, 'LineWidth', 1.5);
            scatter(x, yl(2), 30, col, '^', 'filled', 'MarkerEdgeColor','w');
            text(x, yl(2), sprintf(' %.1f%%', oppo_toeoff_moyen.(cond)), ...
                'Color', col, 'FontSize', 8, 'VerticalAlignment','bottom', 'FontWeight','bold');
        end

        % Opposite heel-strike
        if isfield(oppo_heelstrike_moyen,cond) && ~isnan(oppo_heelstrike_moyen.(cond))
            x = pct2idx(oppo_heelstrike_moyen.(cond), j);
            plot([x x], yl, '-.', 'Color', col, 'LineWidth', 1.5);
            scatter(x, yl(2), 30, col, 's', 'filled', 'MarkerEdgeColor','w');
            text(x, yl(2), sprintf(' %.1f%%', oppo_heelstrike_moyen.(cond)), ...
                'Color', col, 'FontSize', 8, 'VerticalAlignment','middle', 'FontWeight','bold');
        end
    end

    % --- Graduation de 0 à 100% tous les 10% (référence combinée) ---
    xt = round(linspace(1, npts_ref, 11));
    xticks(xt);
    xticklabels(compose('%d%%', 0:10:100));
    xlim([1 npts_ref]);

    grid on;
    if iM == 1 && ~isempty(h_legend)
        legend(h_legend, Conditions, 'Location','best');
    end
end
saveas(gcf, fullfile(output_dir, 'Cycle_Normalise_AllConditions_AllMuscles.png'));
fprintf('Figure globale sauvegardee : Cycle_Normalise_AllConditions_AllMuscles.png\n');

fprintf('\nTermine.\n');