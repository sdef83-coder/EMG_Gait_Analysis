function fft_on_emg_segment(file, EMG_signal, FreqS)
    % Fonction qui applique une FFT sur une plage spécifique du signal EMG
    % Inputs :
    %   - file : le nom du fichier .c3d
    %   - sensor_name : nom du capteur contenant le signal EMG (ex: 'Sensor_3_EMG3')
    %   - FreqS : fréquence d'échantillonnage du signal
    
    % Chargement des données
    data = btkReadAcquisition(file);
    analogs = btkGetAnalogs(data);
    
    % Traitement du signal
    EMGfilt = filtrage(EMG_signal, FreqS, 20, 450);  % Filtrage du signal EMG
    rms_signal = envelop(EMGfilt, FreqS, 0.25); 
    moments = btkGetEvents(data);
    fsrt_frame = btkGetFirstFrame(data);
    lst_frame = btkGetLastFrame(data);
    
    % Tracer le signal EMG filtré et demander à l'utilisateur de sélectionner une plage
    figure;
    plot(EMGfilt);
    title('Sélectionner la plage de données sur le graphique');
    xlabel('Indice');
    ylabel('Amplitude');
    grid on;

    % Utiliser ginput pour sélectionner les points de début et de fin de la plage
    disp('Sélectionnez le point de début et de fin pour la plage de données (cliquez sur deux points sur le graphique)');
    [x, ~] = ginput(2);  % Permet à l'utilisateur de cliquer sur deux points

    start_index = round(x(1));  % Arrondir pour obtenir des indices entiers
    end_index = round(x(2));

    % Extraire la plage de données du signal EMG filtré
    EMGfilt_segment = EMGfilt(start_index:end_index);

    % Appliquer la FFT sur la plage sélectionnée du signal EMG filtré
    N_segment = length(EMGfilt_segment);  % Nombre de points dans le segment
    Y_segment = fft(EMGfilt_segment);  % Transformation de Fourier rapide
    f_segment = (0:N_segment-1)*(FreqS/N_segment);  % Fréquences correspondantes pour la plage
    P2_segment = abs(Y_segment/N_segment);  % Spectre en amplitude
    P1_segment = P2_segment(1:floor(N_segment/2));  % Ne garder que la moitié positive du spectre
    P1_segment(2:end) = 2*P1_segment(2:end);  % Doublez les amplitudes sauf la première valeur

    % Tracer le spectre de fréquence pour la plage sélectionnée
    figure;
    plot(f_segment(1:floor(N_segment/2)), P1_segment);
    title('Spectre de fréquence pour la plage sélectionnée');
    xlabel('Fréquence (Hz)');
    ylabel('Amplitude');
    grid on;
end