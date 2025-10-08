function [EMGenv] = envelope(EMGfilt, FreqEMG, fc, ordre);
    % Fonction pour obtenir l'enveloppe du signal EMG
    %   EMG_raw : signal EMG brut
    %   FreqEMG : fréquence d'échantillonnage du signal EMG
    %   fc : fréquence de coupure pour le filtre passe-bas
    
    % Rectification du signal filtré
    EMGabs = abs(EMGfilt);

    % Normalisation de la fréquence de coupure
    Wn = fc / (FreqEMG / 2);  % Normalisation de la fréquence de coupure
    
    % Conception du filtre passe-bas de second ordre
    [b, a] = butter(ordre, Wn, 'low');
    
    % Application du filtre passe-bas
    EMGenv = filtfilt(b, a, EMGabs);  % Utilisation de la valeur absolue du signal brut
    
end
