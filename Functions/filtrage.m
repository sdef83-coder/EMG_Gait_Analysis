function [EMGfilt] = filtrage(EMG,freq_EMG,fc_low,fc_high)
%filtrage est une fonction qui a pour objectif de filtrer le signal EMG avec un butterworth d'ordre 4 sans apporter de rectification
%   Inputs : EMG (le signal emg brut) ; freq_EMG (la fréquence
%   d'acquisition) ; fc_low (fréquence basse du bandpass) ; fc_high (fréquence haute du bandpass)
Wn = [fc_low,fc_high]/(freq_EMG/2);
[b,a] = butter(4,Wn,"bandpass");
EMGf = filtfilt(b,a,EMG);
EMGfilt = EMGf - mean(EMGf);
end