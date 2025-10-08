clear, close all, clc;

addpath(genpath('C:\Users\silve\OneDrive - Universite de Montreal\Scripts\coherence_analysis'))
cd('C:\Users\silve\OneDrive - Universite de Montreal\Silvere De Freitas - PhD - NeuroBiomech\Scripts\coherence_analysis')

Pairs = {...
    'TAprox','TAdist', 'Sensor_1_EMG1','Sensor_2_EMG2';...
    'GM','SOL','Sensor_4_EMG4','Sensor_3_EMG3';...
    'RF','VL','Sensor_6_EMG6','Sensor_5_EMG5';...
    'VL','ST','Sensor_5_EMG5','Sensor_7_EMG7'} ;

Participants = {'CTL_08','CTL_09'};
Conditions = {'Plat','Medium','High'} ;
Essai = {'01','02','03','04','05','06','07','08','09','10'};

% [cxy,w] = mscohere(x,y,window,noverlap,w) explorer cette fonction
% calcul de coherence entre 2 signaux si ça donne meme chose que
% cross spectrum à discuter avec Fabien

alpha = 0.05 ; % Significant threshold

% Characteristics of the wavelet transformation
FreqMin = 1 ;
FreqMax = 400 ;
Resolution = 1 ;
WaveNumber = 7 ;

tic
for iS = 1%:length(Participants)
    for iC = 1%:length(Conditions)
        for iPair = 1:size(Pairs,1)

            PowSpec_s1 = zeros(400,1000) ; PowSpec_s2 = zeros(400,1000) ; cross_spectrum = zeros(400,1000) ; Ntrials = 0 ;
            % 
            % for iFiles = 1:length(Files) 

                load('SILVERE.mat')
                EMG(:,1) = (c.EMG.Analogs.(Pairs{iPair,3})) ;
                EMG(:,2) = (c.EMG.Analogs.(Pairs{iPair,4})) ;
                
                FreqS = c.c3d.ezc3d.c3d.header.analogs.frameRate ; 
                Cycles = round(c.EMG.Events.Left_Foot_Strike*FreqS)-3000 ; disp('Supprimer le 3000 !!!') 
                Args = WaveletParameters(FreqMin,FreqMax,Resolution,WaveNumber,FreqS) ; % Paramètres de l'ondelette

                [b,a] = butter(2,2*[10 400]/FreqS) ; % EMG cleaning
                EMG = filtfilt(b,a,EMG) ;
                EMG = EMG - repmat(mean(EMG),length(EMG),1) ;

                % calcul carte temps fréquence
                [TFR(:,:,1),period,~,~] = wavelet(EMG(:,1),Args.DT,Args.Pad,Args.DJ,Args.S0,Args.J1,Args.Mother,Args.Cycles) ; % calcul cartes temps-freq
                [TFR(:,:,2),period,~,~] = wavelet(EMG(:,2),Args.DT,Args.Pad,Args.DJ,Args.S0,Args.J1,Args.Mother,Args.Cycles) ;
                Freq = 1./period ;

                Time = linspace(0,length(EMG)/FreqS,length(EMG))' ;                         

                % Découpage cycle par cycle
                for iCycles = 1:length(Cycles)-1
                    TFR_cycle(:,:,1) = TFR(:,Cycles(iCycles):Cycles(iCycles+1),1) ; % découpe les cycles
                    TFR_cycle(:,:,2) = TFR(:,Cycles(iCycles):Cycles(iCycles+1),2) ; % découpe les cycles

                    [X,Y] = meshgrid(1:size(TFR_cycle,2),1:size(TFR_cycle,1)) ; % interpoles les cycles temps-freq sur 1000 points
                    [Xq,Yq] = meshgrid(linspace(1,size(TFR_cycle,2),1000),1:size(TFR_cycle,1)) ;
                    TFR_int(:,:,1) = interp2(X, Y, TFR_cycle(:,:,1), Xq, Yq, 'spline') ;
                    TFR_int(:,:,2) = interp2(X, Y, TFR_cycle(:,:,2), Xq, Yq, 'spline') ;

                    PowSpec_s1 = PowSpec_s1+abs(TFR_int(:,:,1)).^2;
                    PowSpec_s2 = PowSpec_s2+abs(TFR_int(:,:,2)).^2;
                    cross_spectrum = cross_spectrum + (TFR_int(:,:,1)).*conj(TFR_int(:,:,2)) ;

                    Ntrials = Ntrials + 1 ;
                    clear TFR_int TFR_cycle
                end                
        end

            PowSpec_s1 = PowSpec_s1/Ntrials ; % calcul coherence
            PowSpec_s2 = PowSpec_s2/Ntrials ;
            cross_spectrum = abs(cross_spectrum/Ntrials).^2 ;

            Coherence = cross_spectrum./(PowSpec_s1.*PowSpec_s2) ;

            seuil = 1-0.05^(1/(Ntrials-1)) ; % seuil de significativité de la coherence

            Time = linspace(1,1000,100) ;
            figure
            subplot(2,2,1) ; imagesc(Time,Freq,PowSpec_s1) ; title(Pairs{iPair,1})
            subplot(2,2,2) ; imagesc(Time,Freq,PowSpec_s2) ; title(Pairs{iPair,2})
            subplot(2,2,3) ; imagesc(Time,Freq,Coherence>seuil)
            subplot(2,2,4) ; plot(Freq,mean(Coherence,2)) ; hold on % pour les bandes freq calculer l'air sous la courbe de cette fig
            line([0 400],[seuil seuil])
            xlim([1 400])% modifier pour 400hz
            ylim ([0 1])
            %%
            DATA.(['PowSpec_' Pairs{iPair,1}]) = PowSpec_s1 ;
            DATA.(['PowSpec_' Pairs{iPair,2}]) = PowSpec_s2 ;
            DATA.(['Cross_Spectrum_' Pairs{iPair,1} '_' Pairs{iPair,2}]) = cross_spectrum ;
            DATA.(['Coherence_' Pairs{iPair,1} '_' Pairs{iPair,2}]) = Coherence ;
            DATA.Ntrials = Ntrials ;
            DATA.Seuil = seuil ;
            DATA.Freq = Freq ;
          
        end
        save(['C:\Users\fabie\OneDrive - Universite de Montreal\Scripts\Results\' Participants{iP,1} '_' Conditions{iC} '.mat' ],'DATA')
    end
%%
clc, clear, 
load("SILVERE.mat")
CC = c.EMG.Events.Left_Foot_Strike*100
Cycles = round(c.EMG.Events.Left_Foot_Strike*2000)-3000