clear all; close all;

waveFile = 'D:\Down\Kì 5\XLTH\ThucHanh\TH1\TinHieuHuanLuyen\02FVA.wav';
[x, fs]=audioread(waveFile);
frame = x(floor(1*fs) : floor(1.02*fs)); %a
% frame = x(floor(2.3*fs) : floor(2.32*fs)); %e
% frame = x(floor(3.7*fs) : floor(3.72*fs)); %i
% frame = x(floor(5*fs) : floor(5.02*fs)); %o
% frame = x(floor(6.3*fs) : floor(6.32*fs)); %u

formants = test(frame, fs);

function    [A,G,a,r]=autolpc(x,p)
    L=length(x);
    r=[];
    for i=0:p
        r=[r; sum(x(1:L-i).*x(1+i:L))];
    end
    R=toeplitz(r(1:p));
    a=inv(R)*r(2:p+1);
    A=[1; -a];
    G=sqrt(sum(A.*r));
end

function output = myFilter(b, a, x)
    
end

function formants = test(frame, fs)
    frame = frame.*hamming(length(frame));
    preemph=[1 0.63];
    frame = filter(1, preemph, frame);
%     frame = myFilter(1, preemph, frame);
    p = fft(frame);
%     [a, e] = lpc(frame, 12);
    [a,~,~,~] = autolpc(frame,12);
    
    Xlpc = freqz(1, a);
    formant= 20*log10(abs(Xlpc));
    rr = roots(a);
    norm_freq = angle(rr);
    freq_Hz=(norm_freq*fs)/(2*pi);
    freq_hz=sort(abs(freq_Hz));
    v=length(freq_Hz);

    nn = 1;
    for kk = 1:2:length(freq_Hz)
        formants(nn) = freq_hz(kk);
        nn = nn+1;
    end

    disp(formants);
end