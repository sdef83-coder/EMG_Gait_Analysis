%% ANALYSE SPM (STATISTICAL PARAMETRIC MAPPING) - INTER-GROUPES & INTRA-SURFACES
%
% OBJECTIF :
% Ce script réalise deux types d'analyses statistiques sur le cycle de marche :
%   1) INTER-GROUPES : Compare les populations (ex: Adultes vs Adolescents) 
%      sur une même surface (ex: Plat).
%   2) INTRA-GROUPE  : Compare l'effet de la surface (ex: Plat vs High) 
%      au sein d'un même groupe de sujets (test apparié).
%
% ENTRÉES :
%   - Fichiers consolidés : SPM-EMG-ADULTES.mat, SPM-EMG-ADOLESCENTS.mat, etc.
%   - Librairie SPM1D (SnPM pour les tests non-paramétriques).
%
% SORTIES (RÉSULTATS DU SCRIPT) :
%   - IMAGES (.png) : Enregistrées dans /Results/Fig/SPM_EMG/.
%       * Format "SPM_INTER_..." : Différences entre groupes (ex: Adultes vs Ados).
%       * Format "SPM_INTRA_..." : Différences entre surfaces (ex: Plat vs High).
%   - GRAPHIQUES : Doubles tracés (Moyenne+SD en haut / Statistique T en bas).
%   - ZONES SIGNIFICATIVES : Mise en évidence par zones grisées et p-values 
%     affichées directement sur les graphiques pour les zones > seuil critique.
%   - LOG CONSOLE : Suivi en direct des comparaisons significatives trouvées.
%
% -------------------------------------------------------------------------

clc; clear; close all;

%% ===================== CONFIGURATION DES CHEMINS ========================

DATA_DIR = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Matrix';
SPM_PATH = 'C:\Users\silve\Desktop\DOCTORAT\PROGRAMMATION\spm1dmatlab-master';
SAVE_DIR = 'C:\Users\silve\Desktop\DOCTORAT\UNIV MONTREAL\TRAVAUX-THESE\Surfaces Irrégulières\Datas\Script\ActivationMusculaire\Results\Fig\SPM_EMG';

addpath(genpath(SPM_PATH));
if ~exist(SAVE_DIR, 'dir'), mkdir(SAVE_DIR); end
cd(DATA_DIR);

%% ===================== CHARGEMENT DES DONNÉES ===========================

% On charge les 4 groupes. Assurez-vous que les fichiers existent.
GroupData.Adultes       = load('SPM-EMG-ADULTES.mat');
GroupData.Adolescents   = load('SPM-EMG-ADOLESCENTS.mat');
GroupData.Enfants       = load('SPM-EMG-ENFANT.mat');
GroupData.JeunesEnfants = load('SPM-EMG-JEUNESENFANTS.mat');

% Paramètres fixes
GroupNames = fieldnames(GroupData);
Surfaces   = {'Plat', 'Medium', 'High'};
Muscles    = {'EMG_TAprox', 'EMG_TAdist', 'EMG_SOL', 'EMG_GM', 'EMG_VL', 'EMG_RF', 'EMG_ST', 'EMG_GMED'};
alpha      = 0.05;
nPerm      = 5000; % Nombre de permutations pour la robustesse (SnPM)

%% ===================== 1. ANALYSE INTRA-GROUPE (EFFET SURFACE) ==========
% On regarde si, pour les Adolescents (par ex.), marcher sur 'High' change
% l'activation par rapport au 'Plat'.

fprintf('--- DÉBUT : ANALYSE INTRA-GROUPE (PLAT vs SURFACES) ---\n');

% On définit les paires de surfaces à comparer
SurfComps = {'Plat', 'Medium'; 'Plat', 'High'; 'Medium', 'High'};

for iG = 1:length(GroupNames)
    gn = GroupNames{iG}; % Nom du groupe actuel
    
    for iS = 1:size(SurfComps, 1)
        s1 = SurfComps{iS, 1};
        s2 = SurfComps{iS, 2};
        
        for iM = 1:length(Muscles)
            mus = Muscles{iM};
            
            % Extraction (Jambe gauche par défaut)
            Y1 = GroupData.(gn).mean_cycles_combined.(s1).left.(mus);
            Y2 = GroupData.(gn).mean_cycles_combined.(s2).left.(mus);
            
            % Test t Apparié (Paired) : mêmes sujets, deux conditions
            % On ne garde que les sujets ayant des données valides pour les deux surfaces
            valid = ~any(isnan(Y1),2) & ~any(isnan(Y2),2);
            Y1 = Y1(valid,:); Y2 = Y2(valid,:);
            
            if size(Y1,1) < 3, continue; end % Sécurité nb de sujets
            
            % Calcul SPM (Non-paramétrique)
            spm = spm1d.stats.nonparam.ttest_paired(Y1, Y2);
            spmi = spm.inference(alpha, 'two_tailed', true, 'iterations', nPerm);
            
            % Sauvegarde si significatif (p < 0.05)
            if spmi.h0reject
                plot_spm_result(spmi, Y1, Y2, gn, s1, s2, mus, SAVE_DIR, 'INTRA');
            end
        end
    end
end

%% ===================== 2. ANALYSE INTER-GROUPES (EFFET ÂGE) =============
% On regarde si, sur une surface donnée (ex: Medium), les Adolescents 
% sont différents des Adultes.

fprintf('\n--- DÉBUT : ANALYSE INTER-GROUPES (COMPARAISON ÂGES) ---\n');

% Paires de groupes à comparer
GroupComps = {
    'Adultes', 'Adolescents';
    'Adultes', 'Enfants';
    'Adolescents', 'Enfants';
    'Enfants', 'JeunesEnfants'
};

for iS = 1:length(Surfaces)
    surf = Surfaces{iS};
    
    for iC = 1:size(GroupComps, 1)
        g1 = GroupComps{iC, 1};
        g2 = GroupComps{iC, 2};
        
        for iM = 1:length(Muscles)
            mus = Muscles{iM};
            
            % Test t indépendant : deux groupes différents de sujets
            Y1 = GroupData.(g1).mean_cycles_combined.(surf).left.(mus);
            Y2 = GroupData.(g2).mean_cycles_combined.(surf).left.(mus);
            
            % Nettoyage des sujets exclus (NaN)
            Y1(any(isnan(Y1),2),:) = [];
            Y2(any(isnan(Y2),2),:) = [];
            
            if size(Y1,1) < 2 || size(Y2,1) < 2, continue; end
            
            % Calcul SPM (Non-paramétrique)
            spm = spm1d.stats.nonparam.ttest2(Y1, Y2);
            spmi = spm.inference(alpha, 'two_tailed', true, 'iterations', nPerm);
            
            % Sauvegarde si significatif
            if spmi.h0reject
                plot_spm_result(spmi, Y1, Y2, surf, g1, g2, mus, SAVE_DIR, 'INTER');
            end
        end
    end
end

fprintf('\n=== ANALYSE SPM TERMINÉE ===\n');

%% ===================== FONCTION DE VISUALISATION ========================

function plot_spm_result(spmi, Y1, Y2, context, name1, name2, mus, save_dir, type)
    % Prépare une figure invisible pour gagner du temps
    fig = figure('Visible', 'off', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
    
    % En haut : Courbes moyennes +/- SD
    subplot(2,1,1);
    spm1d.plot.plot_meanSD(Y1, 'color', 'b', 'label', name1); hold on;
    spm1d.plot.plot_meanSD(Y2, 'color', 'r', 'label', name2);
    title(sprintf('[%s] %s | %s : %s vs %s', type, context, strrep(mus,'_',' '), name1, name2));
    ylabel('Amplitude EMG Norm.');
    legend('Location','best');
    
    % En bas : Statistique T et seuil de p-value
    subplot(2,1,2);
    spmi.plot();
    spmi.plot_threshold_label();
    spmi.plot_p_values();
    xlabel('% du cycle de marche');
    
    % Sauvegarde PNG
    fname = sprintf('SPM_%s_%s_%s_%s_vs_%s.png', type, context, mus, name1, name2);
    saveas(fig, fullfile(save_dir, fname));
    close(fig);
end