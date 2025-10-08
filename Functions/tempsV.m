function [temps_Vicon] = tempsV(N_c3d,freq_Vicon)
% La fonction Temps permet de dresser un vecteur temps en secondes (échantillon vers secondes)
%   INPUTS : N_Vicon = nombre d'echantillon dans l'essai & freq_Vicon = la fréquence d'acquisition
EMGtempo = (0:N_c3d-1);
temps_Vicon = EMGtempo/freq_Vicon;
end