function [indiceLeft] = indiceLeft(data, analogs, freqS, freqVicon, EMG)
    % La fonction indice va chercher les indices d'events (Heel Strike) de la jambe gauche à partir des mesures Vicon, les met sur une même
    % échelle de temps en seconde que le signal EMG, puis détermine
    % la frame correspondante à chaque Heel Strike de la jambe Gauche.
    %   INPUTS : On rentre data (le lecture btk du fichier c3d), freqS
    %   (frequence d'aquisition de l'EMG), freqVicon (la fréquence
    %   d'acquisition du système Vicon grâce auquel on a les events), EMG (1 signal
    %   électromyographique - peu importe lequel des deux étudiés)
    
    fsrt_frame = btkGetFirstFrame(data);
    moments = btkGetEvents(data);
    N_EMG = length(EMG);
    
    % Récupérer les moments de Heel Strike pour la jambe gauche
    LeftHS = ((moments.Left_Foot_Strike * freqVicon) - fsrt_frame) / freqVicon;
    TempsEMG = tempsE(N_EMG, freqS);  % Temps des échantillons EMG

    % Initialisation des variables
    time_to_find = LeftHS;  
    indices = zeros(size(time_to_find));

    % Recherche des indices les plus proches des moments de Heel Strike dans le signal EMG
    for i = 1:length(time_to_find)
        [~, indices(i)] = min(abs(TempsEMG - time_to_find(i)));
    end

    % Assigner la variable de sortie
    indiceLeft = indices;
end