%% Input signal
[y,Fs] = audioread('D:\Down\Kì 5\XLTH\TinHieuMau\test\LA001.wav'); 

%% Plot wideband and narrowband spec.
subplot(3,1,1);
plot(y);
title('Signal');
%Narrowband : Frame Duration = 30ms , no.fft = 1024 , Frame Shift = 15ms 
%(50% overlap)
subplot(3,1,2);
spec(y,0.03,1024,0.015,Fs);
title('Narrowband');
%Wideband : Frame Duration = 5ms , no.fft = 1024 , Frame Shift = 2.5ms 
%(50% overlap)
subplot(3,1,3);
spec(y,0.005,1024,0.0025 ,Fs);
title('Wideband');

%% Spectrogram display function
function spec(y,WinD,nfft,shiftD,Fs) 
WinL = floor(WinD*Fs); 
% in samples 
L1 = length(y); 
y = y(1:L1); 
% hop size in seconds
shiftL = floor(shiftD*Fs);
nFr = round(length(y)/shiftL); 
%no., of frames
win = hamming(WinL);
 for c = 1:nFr - round(WinL/shiftL) 
    % c : frame count
    FB = (c-1)*shiftL+1; 
    % Beginning Frame 
    FE = FB + WinL -1; 
    %Frame End 
    wseg = y(FB:FE).*win; 
    STFT(:,c) = fft(wseg,nfft);
end
STFTM = abs(STFT(1:nfft/2+1,:)); 
STFTMdb = 20*log10(STFTM);
faxis = (0:nfft/2)*Fs/nfft; 
naxis = (0:size(STFTM,2)-1)*shiftD; 
% in seconds 
 STFTMdbmax = max(STFTMdb(:));
 dbdown = 60;%deciding the range of the plot 
 caxisl = [STFTMdbmax-dbdown STFTMdbmax];% limiting the range of STFT values 
 imagesc(naxis,faxis,STFTMdb,caxisl);
 axis xy; 
 ylabel('Frequency'); 
 xlabel('Length of signal in Seconds'); 
 colorbar;
end