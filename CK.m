clear;
close all;
frameTime = 0.02;
frameShift = 0.01;
F0_Min = 70;
F0_Max = 400;
% path = 'D:\Down\Kì 5\XLTH\CK\TinHieuHuanLuyen\';
% file = { '01MDA', '02FVA', '03MAB', '06FTB' };
path = 'D:\Down\Kì 5\XLTH\CK\TinHieuKiemThu\';
file = { '30FTN', '42FQT', '44MTT', '45MDV' };

i = 1;
% for i = 1 : 4
    figure('Name', char(file(i)));
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    [localSpeech, indexSpeech] = voicedDetection(y, Fs, frameTime, frameShift, path, file(i));
    
    f0 = findF0(y, localSpeech, indexSpeech, Fs, frameTime, frameShift, F0_Max, F0_Min, 2048);
    t = (frameTime - frameShift) : (frameTime - frameShift) : (frameTime - frameShift) * length(f0);
    subplot(5, 1, 3);
    f0 = medianFilter(f0, 3);
    plot(t, f0, '.');
    xlabel('Time(s)');
    ylabel('F0 value(Hz)');
    title('Fundamental frequency (F0)');
    
    [F0_mean, F0_std, deltaMean, deltaStd] = calDelta(f0, path, file(i));
    error = calError(path, file(i), localSpeech);
% end

% Chia frame
% Input:
%   y: mang gia tri cua tin hieu
%   Fs: tan so lay mau (Hz)
%   frameTime: do dai frame (s)
%   frameShift: do lech frame (s)
% Output:
%   data: mang cac frame
%   frameCount: so luong frame
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

% Tinh nang luong ngan han (STE) cua tung frame
% Input:
%   data: mang cac frame
%   frameCount: so luong frame
% Output:
%   ste: mang chua nang luong cua tung frame
function [ste] = shortTimeEnergy(data, frameCount)
    for i = 1:frameCount
        ste(i) = sum(data(i,:).^2);
    end
end

% Tim nguong STE cua tin hieu bang histogram
% Input
%   ste: mang chua nang luong cua tung frame
%   w: tham so tinh nguong, nguoi dung tu khoi tao
% Output:
%   threshold: nguong STE cua tin hieu
function threshold = findThreshold(ste, w)
    [histogram, xValue] = hist(ste, length(ste));
    maxIndex1 = 0;
    maxIndex2 = 0;
    for i = 2 : length(histogram) - 1
        pre = i - 1;
        next = i + 1;
        while (histogram(i) == histogram(next))
            next = next + 1;
        end
        if (histogram(i) > histogram(pre) && histogram(i) > histogram(next))
            if (maxIndex1 == 0)
                maxIndex1 = i;
            else
                maxIndex2 = i;
                break;
            end
        end
    end
    maxHist1 = xValue(maxIndex1);
    maxHist2 = xValue(maxIndex2);
    threshold = (w * maxHist1 + maxHist2) / (w + 1);
end

% Phan biet nguyen am - khoang lang
% Input:
%   y: mang gia tri cua tin hieu
%   Fs: tan so lay mau (Hz)
%   frameTime: do dai frame (s)
%   frameShift: do lech frame (s)
%   path: duong dan chua file tin hieu
%   file: ten file tin hieu
% Output:
%   localSpeech: mang thoi gian thuc cua nguyen am (s)
%   indexSpeech: mang chi so frame cua nguyen am
function [localSpeech, indexSpeech] = voicedDetection(y, Fs, frameTime, frameShift, path, file)
    [data, frameCount] = frameDivide(y, Fs, frameTime, frameShift);
    ste = shortTimeEnergy(data, frameCount);
    ste=ste./max(ste);
    threshold = findThreshold(ste, 10);
    t = (frameTime - frameShift) : (frameTime - frameShift) : (frameTime - frameShift) * length(ste);
    subplot(5, 1, 1);
    plot(t, ste,'r');
    legend('Short time energy');
    title('Short time energy');
    xlabel('Time (s)');
    t = [0 : 1/Fs : length(y)/Fs];
    t = t(1 : end-1);
    subplot(5, 1, 2);
    p1 = plot(t, y);
    hold on;
    [lab, ~, ~] = readFile(path, file);
    x = -1:1;
    for i = 1 : length(lab)
        p2 = plot(lab(i) * ones(size(x)), x, 'r');
    end
    tick = zeros(frameCount, 1);
    for i = 1 : length(ste)
        if (ste(i) >= threshold)
            tick(i) = 1;
        end
    end
    for i = 2 : length(tick) - 1
        if (tick(i) ~= tick(i-1) && tick(i) ~= tick(i+1))
            tick(i) = tick(i-1);
        end
    end
    index = [];
    n = 0;
    for i = 1 : length(tick)
        if (tick(i) == 1)
            n = n + 1;
            index(n) = i;
        end
    end
    % tim diem thoi gian giao nhau giua tieng noi va khoang lang 
    indexSpeech(1) = index(1) - 1; 
    n = 2;
    for i=2:length(index)
        if ((frameTime - frameShift) * (index(i) - index(i-1) + 2) > 0.3)
            indexSpeech(n) = index(i-1);
            indexSpeech(n+1) = index(i) - 1;
            n=n+2;
        end
    end
    indexSpeech(n) = index(i);
    localSpeech = (frameTime - frameShift) * indexSpeech;
    % ve do thi phan doan tieng noi va khoang lang 
    for i = 1:length(localSpeech)
        p3 = plot(localSpeech(i) * ones(size(x)), x,'g');
    end
    hold off; 
    legend([p1, p2, p3], { 'Signal', 'File lab', 'Algorithm' });
    title('Voiced/Unvoiced detection');
    xlabel('Time (s)');
end

% Ham cua so hamming
% Input:
%   M: do dai cua so hamming (mau)
% Output:
%   w: cua so hamming
function w = myHamming(M)
    w = .54 - .46*cos(2*pi*(0:M-1)'/(M-1));
end

% Tim F0 cua tin hieu
% Input:
%   y: mang gia tri cua tin hieu
%   localSpeech: mang thoi gian thuc cua nguyen am (s)
%   indexSpeech: mang chi so frame cua nguyen am
%   Fs: tan so lay mau (Hz)
%   frameTime: do dai frame (s)
%   frameShift: do lech frame (s)
%   F0_Max: tan so gioi han lon nhat (Hz)
%   F0_Min: tan so gioi han nho nhat (Hz)
%   N: so diem tinh fft
% Output:
%   f0: mang cac gia tri f0 (Hz)
function f0 = findF0(y, localSpeech, indexSpeech, Fs, frameTime, frameShift, F0_Max, F0_Min, N)
    frameLength = frameTime * Fs;
    [frames, countF] = frameDivide(y, Fs, frameTime, frameShift);
    p = abs(fft(frames(end - 10, :)', N));
    p = p(1:length(p)/2+1);
    p(2 : end - 1) = 2 * p(2 : end - 1);
    freq = linspace(0, Fs/2, length(p));
    subplot(5, 1, 4);
    plot(freq, p);
    xlabel('Frequency (Hz)');
    ylabel('FFT value');
    legend('Magnitude spectrum');
    title('FFT value of 1 frame unvoice');
    f0 = zeros(countF, 1);
    h = myHamming(frameLength);
    for i = 1 : 2 : length(indexSpeech) - 1
        startIndex = indexSpeech(i) - 1;
        yy = y(floor(localSpeech(i) * Fs) : floor(localSpeech(i+1) * Fs));
        [data, dataCount] = frameDivide(yy, Fs, frameTime, frameShift);
        check = floor(dataCount/2);
        for j = 1 : dataCount
            frame = data(j, :)';
            frame = h.*frame;
            p = abs(fft(frame, N));
            p = p(1:length(p)/2+1);
            p(2 : end - 1) = 2 * p(2 : end - 1);
            freq = linspace(0, Fs/2, length(p));
            if (findpeaks(p, freq, 'NPeaks', 1, 'SortStr', 'descend') > 1)
                [yValue, yPeak] = findpeaks(p, freq, 'MinPeakHeight', 2);
            else
                [yValue, yPeak] = findpeaks(p, freq, 'MinPeakHeight', 0.5);
            end
            if (check == j)
                subplot(5, 1, 5);
                plot(freq, p, yPeak, yValue, 'x');
                xlabel('Frequency (Hz)');
                ylabel('FFT value');
                legend('Magnitude spectrum', 'Harmonics mark');
                title('FFT value of 1 frame voice');
                check = 0;
            end
            tCount = 0;
            tF0 = 0;
            for k = 1 : length(yPeak) - 1
                temp = yPeak(k+1) - yPeak(k);
                if (temp < F0_Max && temp > F0_Min)
                    tCount = tCount + 1;
                    tF0 = tF0 + temp;
                end
            end
            f0(startIndex + j) = tF0 / tCount;
        end
    end
end

% Loc trung vi
% Input:
%   F0: mang cac gia tri f0 (Hz)
%   N: so phan tu lay o hai dau gia tri can tinh trung vi
% Output:
%   F0_med: mang gia tri f0 (Hz) da loc
function F0_med = medianFilter(F0, N)
    for i = 1:N
        F0 = [F0(1); F0; F0(length(F0))];
    end
    F0_med = zeros(length(F0), 1);
    for i = N+1 : length(F0)-N
        t = [];
        for j = -N:N
            t = [t F0(i+j)];
        end
        t = sort(t);
        F0_med(i) = t(round((2*N + 1) / 2));
    end
    F0_med = F0_med(N+1:end-N);
end

% Doc du lieu tu file .lab
% Input:
%   path: duong dan chua file tin hieu
%   file: ten file tin hieu
% Output:
%   lab: mang chua cac gia tri thoi gian cua nguyen am (s)
%   meanLab: gia tri F0_mean cua file .lab (Hz)
%   stdLab: gia tri F0_std cua file .lab (Hz)
function [lab, meanLab, stdLab] = readFile(path, file)
    filename = char(strcat(path, file, '.lab'));
    f = fopen(filename, 'r');
    line = {};
    i = 0;
    while ~feof(f)
       i = i + 1;
       line(i) = cellstr(fgets(f));
    end
    k = 1;
    for j = 2 : 2 : 10
        vowel = split(line{j});
        lab(k) = str2double(vowel(1));
        lab(k+1) = str2double(vowel(2));
        k = k + 2;
    end
    meanLab = split(line(i-1));
    meanLab = str2double(meanLab(2));
    stdLab = split(line(i));
    stdLab = str2double(stdLab(2));
    fclose(f);
end

% Tinh sai so cua thuat toan so voi file .lab
% Input:
%   F0: mang cac gia tri f0 (Hz)
%   path: duong dan chua file tin hieu
%   file: ten file tin hieu
% Output:
%   F0_mean: gia tri F0 trung binh cua thuat toan (Hz)
%   F0_std: gia tri do lech chuan F0 cua thuat toan (Hz)
%   deltaMean: gia tri do lech F0 trung binh cua thuat toan (Hz)
%   deltaStd: gia tri do lech do lech chuan F0  cua thuat toan (Hz)
function [F0_mean, F0_std, deltaMean, deltaStd] = calDelta(F0, path, file)
    F0_mean = 0;
    count = 0;
    for i = 1 : length(F0)
        if F0(i) > 0
            count = count + 1;
            F0_mean = F0_mean + F0(i);
        end
    end
    F0_mean = round(F0_mean / count, 2);
    F0_std = 0;
    for i = 1 : length(F0)
        if F0(i) > 0
            F0_std = F0_std + (F0(i) - F0_mean)^2;
        end
    end
    F0_std = round(sqrt(F0_std / count), 2);
    [~, meanLab, stdLab] = readFile(path, file);
    deltaMean = abs(F0_mean - meanLab);
    deltaStd = abs(F0_std - stdLab);
end

% Tinh sai so tuyet doi cua bien thoi gian (ms)
% Input:
%   path: duong dan chua file tin hieu
%   file: ten file tin hieu
%   localSpeech: mang thoi gian thuc cua nguyen am (s)
% Output:
%   error: mang gia tri sai so
function error = calError(path, file, localSpeech)
    [lab, ~, ~] = readFile(path, file);
    error = [];
    for i = 1 : length(lab)
        error(i) = abs(lab(i) - localSpeech(i));
    end
    error = error.*1000;
end





