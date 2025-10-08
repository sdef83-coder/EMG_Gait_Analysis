function [temps_EMG] = tempsE(N_EMG,freq_EMG)
% La fonction Temps permet de dresser un vecteur temps en secondes (échantillon vers secondes)
%   INPUTS : N_EMG = nombre d'echantillon dans l'essai & freq_EMG = la fréquence d'acquisition
EMGtempo = (0:N_EMG-1);
temps_EMG = EMGtempo/freq_EMG;
end