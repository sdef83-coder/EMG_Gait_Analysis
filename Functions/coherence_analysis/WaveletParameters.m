function [Args] = WaveletParameters(FreqMin,FreqMax,Resolution,WaveNumber,FreqS)
% Define arguments to compute TFR from FReqMin:Resolution:FreqMax
% WaveNumber: the number of cycles
% FreqS: Frequency sampling

fourier_factor = single((4*pi)/(WaveNumber + sqrt(2 + WaveNumber^2))) ; 
S0 = FreqMin*fourier_factor ;
J1 = FreqMax*fourier_factor ;
nvoice = length(FreqMin:Resolution:FreqMax)/(J1-S0+1) ;

% Arguments of the wavelet transform
Args=struct('DT',1/FreqS,...
    'Pad',1,...      % pad the time series with zeroes (recommended)
    'S0',S0,...   smallest scale of the wavelet.  Default is 2*DT.
    'J1',J1,...  the # of scales minus one S0 up to S0*2^(J1*DJ)
    'DJ',nvoice, ...    spacing between discrete scales
    'Mother','Morlet',...
    'MaxScale',[],...
    'Cycles',WaveNumber) ;