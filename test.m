clear;
close all;
frameTime = 0.02;
frameShift = 0.01;
path = 'D:\Down\Kì 5\XLTH\ThucHanh\TH1\TinHieuHuanLuyen\';
file = { '01MDA', '02FVA', '03MAB', '06FTB', '30FTN', '42FQT', '44MTT', '45MDV' };



% [~, freq, time, power] = spectrogram(y, floor(5*10^(-3)*Fs), floor(2*10^(-3)*Fs), 2048, Fs, 'yaxis');
% [yP, xP] = findpeaks(power(:,262)', freq');

i = 8;
k = 101;
% for i = 1 : length(file)
    figure('Name', char(file(i)));
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    
    [data, frameCount] = frameDivide(y, Fs, frameTime, frameShift);
    
    h = hamming(frameTime * Fs);
    frame = data(k, :)';
    frame = h.*frame;
%     p = abs(fft(frame));
%     freq = linspace(0, Fs/2, length(p));
%     p = log(p);
%     s = abs(ifft(p));
%     subplot(3,1,1);
%     plot(freq, s);
    C = real(ifft(log(abs(fft(frame)))));
    C = C(1:length(C)/2+1);
    subplot(3,1,1);
    freq = linspace(0, Fs/2, length(C));
%     plot(freq, C);
    plot(C);
    s1 = ifft(C);
    s1 = fft(s1);
    subplot(3,1,2);
    plot(abs(s1));
% end

function [data, frameCount] = frameDivide(y, Fs, frameTime, frameShift)
    N = length(y);
    frameLength = frameTime * Fs;
    shiftLength = frameShift * Fs;
    frameCount = floor((N - frameShift * Fs) / (frameLength - shiftLength));
    temp = 0;
    for i = 1 : frameCount
        data(i, :) = y(temp + 1 : temp + frameLength);
        temp = temp + shiftLength;
    end
end

