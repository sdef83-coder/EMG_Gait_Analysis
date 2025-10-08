function [signal_cleaned, EMGenvcleaned] = cleanEMGdata(EMGfilt, FreqS, method, thresholdFactor, minDuration, visualize)
    % Fonction pour nettoyer le signal EMG en détectant et traitant les artefacts
    % Entrées :
    %   EMGfilt - Signal EMG filtré
    %   FreqS - Fréquence d'échantillonnage
    %   method - Méthode pour détecter les artefacts ('STD', 'MAD', 'IQR', 'CUSTOM')
    %   thresholdFactor - Facteur de seuil pour la détection des artefacts
    %   minDuration - Durée minimale des artefacts (en échantillons)
    %   visualize - Booléen pour afficher les graphiques (true/false)
    %
    % Sorties :
    %   signal_cleaned - Signal nettoyé (artefacts supprimés et interpolés)
    %   EMGenvcleaned - Enveloppe du signal nettoyé sans pics isolés

    % 1. Détection des artefacts
    artifacts_info = characterizeArtifacts(EMGfilt, FreqS, 'Method', method, 'ThresholdFactor', thresholdFactor, 'MinDuration', minDuration, 'Visualize', visualize);
    
    % Obtention des seuils calculés
    upper_threshold = artifacts_info.upper_threshold;
    lower_threshold = artifacts_info.lower_threshold;
    
    % Détection des valeurs aberrantes (outliers)
    outliers = (EMGfilt > upper_threshold) | (EMGfilt < lower_threshold);
    
    % 2. REMPLACEMENT DES ARTEFACTS
    % Remplacer les artefacts par des NaN pour interpolation ultérieure
    signal_cleaned = EMGfilt;
    signal_cleaned(outliers) = NaN;
    
    % Interpolation pour remplacer les NaN (artefacts)
    signal_cleaned = fillmissing(signal_cleaned, 'pchip');
    
    % 3. ÉLIMINATION DES PICS ISOLÉS - Filtre médian (pour les pics très courts)
    window_size = 5; % Taille de fenêtre pour le filtre médian (ajustable)
    signal_cleaned = medfilt1(signal_cleaned, window_size);
    
    % 4. Calcul de l'enveloppe des signaux filtrés
    abscleaned = abs(signal_cleaned);
    EMGenvcleaned = envelope(abscleaned, FreqS, 9, 2);
    