function [rms_signal] = envelop(EMGfilt,freq_EMG,window)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
EMGabs = abs(EMGfilt);
fenetre = round(window * freq_EMG); % Calcul de la taille de la fenêtre
rms_signal = sqrt(movmean(EMGabs.^2, fenetre)); %Calcul du RMS glissant
end