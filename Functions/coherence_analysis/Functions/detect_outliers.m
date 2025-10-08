function [outlier_indices] = detect_outliers(muscle_maxima, cycle_offset, threshold)
    % Cette fonction détecte les outliers à partir des maxima de cycles et les seuils basés sur les écarts-types.
    % Paramètres d'entrée :
    %   - muscle_maxima : Liste des maxima des cycles pour un muscle donné
    %   - cycle_offset : Décalage pour la concaténation des cycles (pour le tracé)

    % Calcul de la moyenne et de l'écart-type des maxima des cycles
    mean_max = mean(muscle_maxima);
    std_max = std(muscle_maxima);

    % Détecter les outliers (cycles dont la valeur est supérieure à 2 écarts-types de la moyenne)
    outlier_upper_threshold = mean_max + threshold*std_max;  % Seulement les cycles dont les valeurs sont supérieures à 2 écarts-types sont considérés comme outliers

    % Identification des outliers
    outlier_indices = (muscle_maxima > outlier_upper_threshold);  % 1 si outlier, 0 sinon

    % Affichage des seuils
    plot([0, cycle_offset], [outlier_upper_threshold, outlier_upper_threshold], 'r--', 'LineWidth', 1.5);  % Ligne du seuil supérieur
end