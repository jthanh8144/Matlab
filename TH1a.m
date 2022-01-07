clear;
close all;
frameTime = 0.02;
frameShift = 0.01;
thresholdACF = 0.3;
F0_Min = 70;
F0_Max = 400;
% path = 'D:\Down\Kì 5\XLTH\ThucHanh\TH1\TinHieuHuanLuyen\';
% file = { '01MDA', '02FVA', '03MAB', '06FTB', '30FTN', '42FQT', '44MTT', '45MDV' };
path = 'D:\Down\Kì 5\XLTH\ThucHanh\TH1\TinHieuKiemThu-44k\';
file = { '04MHB', '05MVB', '07FTC', '08MLD', '09MPD', '10MSD', '12FTD', '14FHH', '16FTH', '24FTL' };

resultACF = {};
resultFFT = {};
F0_meanACF = zeros;
F0_stdACF = zeros;
F0_meanFFT = zeros;
F0_stdFFT = zeros;
% i = 2;
for i = 1 : length(file)
    figure('Name', char(file(i)));
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    [localSpeech, indexSpeech] = voicedDetection(y, Fs, frameTime, frameShift, file(i));
    
    f0ACF = findF0usingACF(y, indexSpeech, Fs, frameTime, frameShift, F0_Min, F0_Max, thresholdACF);
    t = (frameTime - frameShift) : (frameTime - frameShift) : (frameTime - frameShift) * length(f0ACF);
    [F0_mean, F0_std, ~, ~] = calDelta(f0ACF, path, file(i));
    f0ACF = filterF0(f0ACF, F0_mean, F0_std);
    subplot(5, 2, 3:4);
    f0ACF = medianFilter(f0ACF, 2);
    plot(t, f0ACF, '.');
    xlabel('Time(s)');
    ylabel('F0 value(Hz)');
    title('Fundamental frequency (F0) with ACF');
%     [F0_meanACF(i), F0_stdACF(i), deltaMean, deltaStd] = calDelta(f0ACF, path, file(i));
    [F0_meanACF(i), F0_stdACF(i), ~, ~] = calDelta(f0ACF, path, file(i));
    resultACF(i) = cellstr(genderDetection(F0_meanACF(i)));
    
    f0 = findF0usingFFT(y, indexSpeech, Fs, frameTime, frameShift, F0_Max, F0_Min, 2048);
    t = (frameTime - frameShift) : (frameTime - frameShift) : (frameTime - frameShift) * length(f0);
    [F0_mean, F0_std, ~, ~] = calDelta(f0, path, file(i));
    f0 = filterF0(f0, F0_mean, F0_std);
    subplot(5, 2, 5:6);
    f0 = medianFilter(f0, 2);
    plot(t, f0, '.');
    xlabel('Time(s)');
    ylabel('F0 value(Hz)');
    title('Fundamental frequency (F0) with FFT');
    
%     [F0_meanFFT(i), F0_stdFFT(i), deltaMean, deltaStd] = calDelta(f0, path, file(i));
    [F0_meanFFT(i), F0_stdFFT(i), ~, ~] = calDelta(f0, path, file(i));
    resultFFT(i) = cellstr(genderDetection(F0_meanFFT(i)));
end
clc;
disp(file);
disp(F0_meanACF);
% disp(F0_stdACF);
disp(resultACF);
% disp(F0_meanFFT);
% disp(F0_stdFFT);
% disp(resultFFT);
disp(accuracyGenderDetection(file, resultACF));
% disp(accuracyGenderDetection(file, resultFFT));

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
function threshold = findThresholdSTE(ste, w)
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
%   file: ten file tin hieu
% Output:
%   localSpeech: mang thoi gian thuc cua nguyen am (s)
%   indexSpeech: mang chi so frame cua nguyen am
function [localSpeech, indexSpeech] = voicedDetection(y, Fs, frameTime, frameShift, file)
    [data, frameCount] = frameDivide(y, Fs, frameTime, frameShift);
    ste = shortTimeEnergy(data, frameCount);
    ste=ste./max(ste);
    threshold = findThresholdSTE(ste, 10);
    t = [0 : 1/Fs : length(y)/Fs];
    t = t(1 : end-1);
    subplot(5, 2, 1:2);
    plot(t, y);
    xlabel('Time(s)');
    ylabel('Signal');
    title(['Input signal: ', char(file)]);
    
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
end

% Ham tu tuong quan
% Input:
%   y: mang gia tri cua tin hieu
%   frame_length: Do dai mot khung (mau)
% Output:
%   yy: mang gia tri tu tuong quan
function yy = ACFunc(y, frame_length)
    yy = zeros(1, round(frame_length));
    for i = 0:length(y)-1
        yy(i+1)= 0;
        for j = 1:length(y)-i
            yy(i+1) = yy(i+1) + y(j)*y(i+j);
        end
    end
end

% Tim f0 cua mot khung bang ACF
% Input:
%   yy: mang tu tuong quan
%   Fs: tan so lay mau (Hz)
%   F0_Min: tan so gioi han nho nhat
%   F0_Max: tan so gioi han lon nhat
%   threshold: nguong
% Output:
%   f0: gia tri f0 cua khung
function f0 = find_F0_1frame(yy, Fs, F0_Min, F0_Max, threshold)
    t_Min = Fs / F0_Max;
    t_Max = Fs / F0_Min;
    f0 = 0;
    maxA = [];
    for i = floor(t_Min):floor(t_Max)
        if(yy(i+1) < yy(i) && yy(i-1) < yy(i) && yy(i) > threshold*yy(1))
            maxA = [maxA i];
        end
    end
    t = 0;
    for i=1:length(maxA)
        if (t == 0) 
            t = maxA(i);
        elseif (yy(maxA(i)) > yy(t))
            t = maxA(i);
        end
    end
    if (t == 0)
        f0 = 0;
    else
        f0 = Fs / (t - 1);
    end
end

% Tim F0 cua tin hieu bang ACF
% Input:
%   y: mang gia tri cua tin hieu
%   Fs: tan so lay mau (Hz)
%   frame_time: do dai cua mot khung (s)
%   F0_Min: tan so gioi han nho nhat
%   F0_Max: tan so gioi han lon nhat
%   threshold: nguong
% Output:
%   F0: mang gia tri F0 (Hz)
function f0 = findF0usingACF(y, indexSpeech, Fs, frameTime, frameShift, F0_Min, F0_Max, threshold)
    frameLength = frameTime * Fs;
    [frames, countF] = frameDivide(y, Fs, frameTime, frameShift);
%     Ve ket qua trung gian
    yy = ACFunc(frames(end,:)', frameLength);
    n=1:length(yy);
    t = n/Fs;
    subplot(5,2,7);
    plot(t, yy);
    xlabel('Time(s)');
    ylabel('ACF value');
    title('ACF value of 1 frame unvoice');
    
    f0 = zeros(countF, 1);
    for i = 1 : 2 : length(indexSpeech) - 1
        check = indexSpeech(i) + 10;
        for j = indexSpeech(i) : indexSpeech(i + 1)
            yt = frames(j, :)';
            yy = ACFunc(yt, frameLength);
            f0(j) = find_F0_1frame(yy, Fs, F0_Min, F0_Max, threshold);
%             Ve ket qua trung gian
            if (check == j)
                n=1:length(yy);
                t = n/Fs;
                subplot(5,2,8);
                plot(t, yy);
                xlabel('Time(s)');
                ylabel('ACF value');
                title('ACF value of 1 frame voice');
            end
        end
    end
end

% Ham cua so hamming
% Input:
%   M: do dai cua so hamming (mau)
% Output:
%   w: cua so hamming
function w = myHamming(M)
    w = .54 - .46*cos(2*pi*(0:M-1)'/(M-1));
end

% Tim F0 cua tin hieu bang FFT
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
function f0 = findF0usingFFT(y, indexSpeech, Fs, frameTime, frameShift, F0_Max, F0_Min, N)
    frameLength = frameTime * Fs;
    [frames, countF] = frameDivide(y, Fs, frameTime, frameShift);
    f0 = zeros(countF, 1);
    h = myHamming(frameLength);
    
%     Ve ket qua trung gian
    p = abs(fft(h.*frames(end - 10, :)', N));
    p = p(1:length(p)/2+1);
    p(2 : end - 1) = 2 * p(2 : end - 1);
    freq = linspace(0, Fs/2, length(p));
    subplot(5, 2, 9);
    plot(freq, p);
    xlabel('Frequency (Hz)');
    ylabel('FFT value');
    legend('Magnitude spectrum');
    title('FFT value of 1 frame unvoice');
    
    for i = 1 : 2 : length(indexSpeech) - 1
        check = indexSpeech(1) + 10;
        for j = indexSpeech(i) : indexSpeech(i + 1)
            frame = frames(j, :)';
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
%             Ve ket qua trung gian
            if (check == j)
                subplot(5, 2, 10);
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
            f0(j) = tF0 / tCount;
        end
    end
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
    deltaMean = 0;
    deltaStd = 0;
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
%     [~, meanLab, stdLab] = readFile(path, file);
%     deltaMean = abs(F0_mean - meanLab);
%     deltaStd = abs(F0_std - stdLab);
end

% Loc cac F0 qua sai
% Input:
%   F0: mang cac gia tri f0 (Hz)
%   F0_mean: gia tri F0 trung binh (Hz)
%   F0_std: gia tri do lech chuan F0 (Hz)
% Ouput:
%   F0_filter: mang gia tri f0 (Hz) da loc
function F0_filter = filterF0(F0, F0_mean, F0_std)
    F0_filter = zeros(length(F0), 1);
    for i = 1 : length(F0)
        if (F0(i) > 0)
            if (abs(F0(i) - F0_mean) > F0_std) 
                F0_filter(i) = F0_mean;
            else
                F0_filter(i) = F0(i);
            end
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

% Phan biet gioi tinh giong noi bang F0
% Input:
%   F0_mean: gia tri F0 trung binh (Hz)
% Output:
%   result: gioi tinh cua tin hieu
function result = genderDetection(F0_mean)
    if (F0_mean >= 70 && F0_mean <= 150)
        result = 'Male';
    else
        result = 'Female';
    end
end

% Tinh do chinh xac phan biet gioi tinh
% Input:
%   file: mang ten file tin hieu
%   gender: mang gioi tinh phan biet duoc
% Output:
%   result: do chinh xac phan biet gioi tinh (%)
function result = accuracyGenderDetection(file, gender)
    reference = {};
    for i = 1 : length(file)
        if (file{i}(3) == 'M')
            reference(i) = cellstr('Male');
        else
            reference(i) = cellstr('Female');
        end
    end
    count = 0;
    for i = 1 : length(gender)
        if (strcmp(reference(i), gender(i)) == true)
            count = count + 1;
        end
    end
    result = count / length(gender) * 100;
end









