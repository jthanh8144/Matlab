%
% pr9_2_1
clear all; close all;

% waveFile='snn27.wav';
waveFile = 'D:\Down\Kì 5\XLTH\ThucHanh\TH1\TinHieuHuanLuyen\45MDV.wav';
[x, fs]=audioread(waveFile);
x = x(floor(1.03*fs) : floor(1.06*fs));
t = 0 : 1/fs: length(x)/fs - 1/fs;
subplot(2,1,1);
plot(t, x);
u=filter([1 -.99],1,x);                          
wlen=length(u);                                 
cepstL=6;                                        
wlen2=wlen/2;               
freq=(0:wlen2-1)*fs/wlen;                      
u2=u.*hamming(wlen);		                     
U=fft(u2);                                       
U_abs=log(abs(U(1:wlen2)));                      
Cepst=ifft(U_abs);                               
cepst=zeros(1,wlen2);           
cepst(1:cepstL)=Cepst(1:cepstL);                 
cepst(end-cepstL+2:end)=Cepst(end-cepstL+2:end);
spect=real(fft(cepst));                          

%[Loc,Val]=findpeaks(spect);                      
[Val,Loc]=findpeaks(spect);
FRMNT=freq(Loc);                                 
% ×÷Í¼
subplot(2,1,2);
pos = get(gcf,'Position');
set(gcf,'Position',[pos(1), pos(2)-100,pos(3),(pos(4)-140)]);
plot(freq,U_abs,'k'); 
hold on; axis([0 4000 -6 2]); grid;
plot(freq,spect,'k','linewidth',2); 
xlabel('ÆµÂÊ/Hz'); ylabel('·ùÖµ/dB');
title('ÐÅºÅÆµÆ×,°üÂçÏßºÍ¹²Õñ·åÖµ')
% fprintf('%5.2f   %5.2f   %5.2f   %5.2f\n',FRMNT);
disp(FRMNT);
for k=1 : 4
    plot(freq(Loc(k)),Val(k),'kO','linewidth',2);
    line([freq(Loc(k)) freq(Loc(k))],[-6 Val(k)],'color','k',...
        'linestyle','-.','linewidth',2);
end
