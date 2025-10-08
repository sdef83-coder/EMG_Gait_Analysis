%% COMPARAISON SPM EMG ENFANTS VS JEUNES ENFANTS VS ADULTES

clc; clear; close all;

% Chemins
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Matrix');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\PROGRAMMATION\spm1dmatlab-master'));

% Chargement des données
GroupData.Enfants = load('SPM-EMG-Enfant.mat');
GroupData.Adultes = load('SPM-EMG-ADULTES.mat');
GroupData.JeunesEnfants = load('SPM-EMG-JeunesEnfants.mat');

% Comparaisons à effectuer
Comparisons = {
    'Adultes', 'Enfants';
    'Adultes', 'JeunesEnfants';
    'JeunesEnfants', 'Enfants'
};

% Paramètres
saveFolder = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Fig\SPM_EMG';
if ~exist(saveFolder, 'dir'); mkdir(saveFolder); end

x = 1:100;
muscles = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'};
conditions = {'Plat', 'Medium', 'High'};
alpha = 0.05;
affiche_avec_std = @(mean_data, std_data, color) fill([x fliplr(x)], [mean_data + std_data, fliplr(mean_data - std_data)], color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

% Boucle sur chaque comparaison
for comp = 1:size(Comparisons, 1)
    group1 = Comparisons{comp, 1};
    group2 = Comparisons{comp, 2};

    for m = 1:length(muscles)
        muscle = muscles{m};
        fig = figure('Position',[100 100 1400 800]);
        sgtitle(['Comparaison ' group1 ' vs ' group2 ' - EMG (SnPM) - Muscle : ' muscle]);

        for c = 1:length(conditions)
            cond = conditions{c};

            try
                Y1 = GroupData.(group1).mean_cycles_combined.(cond).(muscle);
                Y2 = GroupData.(group2).mean_cycles_combined.(cond).(muscle);
            catch
                warning(['❌ Données manquantes pour ' muscle ' - ' cond]);
                continue
            end

            if size(Y1,2) ~= 100 || size(Y2,2) ~= 100
                warning(['❌ Format incorrect pour ' muscle ' - ' cond]);
                continue
            end

            % Moyenne et écart-type
            mean1 = mean(Y1, 1);
            std1 = std(Y1, 0, 1);
            mean2 = mean(Y2, 1);
            std2 = std(Y2, 0, 1);

            % Tracé EMG
            subplot(length(conditions), 2, (c-1)*2+1);
            hold on;
            affiche_avec_std(mean1, std1, 'b');
            affiche_avec_std(mean2, std2, 'r');
            plot(x, mean1, 'b', 'LineWidth', 1.5);
            plot(x, mean2, 'r', 'LineWidth', 1.5);
            title(['EMG - ' cond]);
            legend(group1, group2);
            xlabel('% cycle de marche');
            ylabel('Activité EMG (norm.)');
            grid on;

            % Analyse SnPM
            subplot(length(conditions), 2, (c-1)*2+2);
            hold on;
            disp(['SnPM pour ' muscle ' - ' cond ' : ' group1 ' vs ' group2]);
            spm = spm1d.stats.nonparam.ttest2(Y1, Y2);
            nPermMax = spm.permuter.nPermTotal;  % nombre max possible de permutations
            iterations = min(5000, nPermMax);   % 5000 si possible, sinon le max
            spmi = spm.inference(alpha, 'two_tailed', true, 'interp', true, 'iterations', iterations);
            spmi.plot();
            spmi.plot_threshold_label();
            spmi.plot_p_values();
            title(['SnPM1D - ' cond]);
            ylabel('T-statistique');
            xlabel('% cycle de marche');
            grid on;

            % Clusters significatifs
            if spmi.h0reject && isstruct(spmi.clusters)
                for i = 1:length(spmi.clusters)
                    if spmi.clusters(i).h0reject
                        cluster = spmi.clusters(i).indices;
                        fill([x(cluster) fliplr(x(cluster))], ...
                             [max(ylim)*ones(size(cluster)) fliplr(min(ylim)*ones(size(cluster)))], ...
                             'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                        text(mean(x(cluster)), max(ylim)*0.9, '*', ...
                             'FontSize', 16, 'Color', 'k', 'FontWeight', 'bold', ...
                             'HorizontalAlignment', 'center');
                    end
                end
            else
                disp(['→ Aucune différence significative pour ' muscle ' - ' cond]);
            end
        end

        % Sauvegarde figure
        saveName = fullfile(saveFolder, ['SnPM_EMG_' muscle '_' group1 '_vs_' group2]);
        exportgraphics(fig, [saveName '.png'], 'Resolution', 300);
        close(fig);
    end
end

%% COMPARAISON INTRAGROUPE - Plat vs High - EMG SnPM
clc; clear; close all;

% Chemins
cd('C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Matrix');
addpath(genpath('C:\Users\silve\Desktop\DOCTORAT\PROGRAMMATION\spm1dmatlab-master'));

Enfants = load('SPM-EMG-Enfant.mat');
Adultes = load('SPM-EMG-Adultes.mat');
Jeunes_Enfants = load('SPM-EMG-JeunesEnfants.mat');

saveFolder = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Fig\SPM_EMG_CONDITIONS';
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder)
end

x = 1:100;
muscles = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'};
alpha = 0.05;
groupes = {'Adultes', 'Enfants', 'Jeunes_Enfants'};

% Fonction affichage bande ±1 SD
affiche_avec_std = @(mean_data, std_data, color) fill([x fliplr(x)], ...
    [mean_data + std_data, fliplr(mean_data - std_data)], ...
    color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');

for g = 1:length(groupes)
    nomGroupe = groupes{g};
    switch nomGroupe
    case 'Adultes'
        DATA = Adultes;
    case 'Enfants'
        DATA = Enfants;
    case 'Jeunes_Enfants'
        DATA = Jeunes_Enfants;
    otherwise
        error(['Groupe inconnu : ' nomGroupe]);
end

    for m = 1:length(muscles)
        muscle = muscles{m};

        fig = figure('Position',[100 100 1400 600]);
        sgtitle([nomGroupe ' - Comparaison EMG Plat vs High - Muscle : ' muscle]);

        % Extraction
        try
            Y1 = DATA.mean_cycles_combined.Plat.(muscle);   % Plat
            Y2 = DATA.mean_cycles_combined.High.(muscle);   % High
        catch
            warning(['❌ Données manquantes pour ' nomGroupe ' - ' muscle]);
            close(fig);
            continue;
        end

        if size(Y1,2) ~= 100 || size(Y2,2) ~= 100
            warning(['❌ Format incorrect pour ' nomGroupe ' - ' muscle]);
            close(fig);
            continue;
        end

        % EMG signal (± SD)
        subplot(1,2,1);
        hold on
        mean1 = mean(Y1, 1); std1 = std(Y1, 0, 1);
        mean2 = mean(Y2, 1); std2 = std(Y2, 0, 1);

        affiche_avec_std(mean1, std1, 'b');
        affiche_avec_std(mean2, std2, 'r');
        plot(x, mean1, 'b', 'LineWidth', 1.5);
        plot(x, mean2, 'r', 'LineWidth', 1.5);
        title(['EMG - Plat vs High']);
        legend('Plat', 'High');
        xlabel('% cycle de marche');
        ylabel('Activité EMG (norm.)');
        grid on;

        % SnPM
        subplot(1,2,2);
        hold on
        disp(['SnPM - ' nomGroupe ' - ' muscle ' : Plat vs High']);
        spm = spm1d.stats.nonparam.ttest2(Y1, Y2);
        nPermMax = spm.permuter.nPermTotal;  % nombre max possible de permutations
        iterations = min(5000, nPermMax);   % 5000 si possible, sinon le max
        spmi = spm.inference(alpha, 'two_tailed', true, 'interp', true, 'iterations', iterations);
        spmi.plot();
        spmi.plot_threshold_label();
        spmi.plot_p_values();
        title('SnPM1D - Plat vs High');
        ylabel('T-statistique');
        xlabel('% cycle de marche');
        grid on;

        % Clusters
        if spmi.h0reject && isstruct(spmi.clusters)
            for i = 1:length(spmi.clusters)
                if spmi.clusters(i).h0reject
                    cluster = spmi.clusters(i).indices;
                    fill([x(cluster) fliplr(x(cluster))], ...
                        [max(ylim)*ones(size(cluster)) fliplr(min(ylim)*ones(size(cluster)))], ...
                        'k', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
                    text(mean(x(cluster)), max(ylim)*0.9, '*', ...
                         'FontSize', 16, 'Color', 'k', 'FontWeight', 'bold', ...
                         'HorizontalAlignment', 'center');
                end
            end
        else
            disp(['→ Aucune différence significative pour ' nomGroupe ' - ' muscle]);
        end

        % Sauvegarde
        saveName = fullfile(saveFolder, ['SnPM_' nomGroupe '_Plat_vs_High_' muscle]);
        exportgraphics(fig, [saveName '.png'], 'Resolution', 300);
        close(fig);
    end
end