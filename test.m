clear;
close all;
frameTime = 0.02;
frameShift = 0.01;
N = 13; % 13 26 39
K = 4; % 2 -> 5
NFFT = 1024;
path = 'D:\Down\K? 5\XLTH\ThucHanh\TH1\NguyenAmHuanLuyen-16k\';
folder = { '01MDA', '02FVA', '03MAB', '04MHB', '05MVB', '06FTB', '07FTC', '08MLD', '09MPD', '10MSD', '11MVD', '12FTD', '14FHH', ...
                '15MMH', '16FTH', '17MTH', '18MNK', '19MXK', '20MVK', '21MTL', '22MHL' };
pathKT = 'D:\Down\K? 5\XLTH\ThucHanh\TH1\NguyenAmKiemThu-16k\';
folderKT = { '23MTL', '24FTL', '25MLM', '27MCM', '28MVN', '29MHN', '30FTN', '32MTP', '33MHP', '34MQP', '35MMQ', '36MAQ', ...
                    '37MDS', '38MDS', '39MTS', '40MHS', '41MVS', '42FQT', '43MNT', '44MTT', '45MDV' };
file = { 'a', 'e', 'i', 'o', 'u' };

% uncomment line 18 when training, when use line 18 comment line 21-22
% [fingerprints, fftDB] = getFingerprints(frameTime, frameShift, N, K, path, folder, file, NFFT);
getFingerprints(frameTime, frameShift, N, K, path, folder, file, NFFT);

% uncomment line 21 -> 24 when testing
fingerprints = readmatrix(char(strcat(path, 'fingerprints.xlsx')));
fftDB = readmatrix(char(strcat(path, 'fftDB.xlsx')));
[result, accuracy, resultFFT, accuracyFFT] = runTesting(fingerprints, fftDB, pathKT, folderKT, file, frameTime, frameShift, N, K, NFFT);
plotFigure(fingerprints, fftDB, K, file);


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
        temp = temp + (frameLength - shiftLength);
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
function [localSpeech, indexSpeech] = voicedDetection(y, Fs, frameTime, frameShift)
    [data, frameCount] = frameDivide(y, Fs, frameTime, frameShift);
    ste = shortTimeEnergy(data, frameCount);
    ste=ste./max(ste);
    threshold = findThresholdSTE(ste, 10);
    
    tick = zeros(frameCount, 1);
    for i = 1 : length(ste)
        if (ste(i) >= threshold)
            tick(i) = 1;
        end
    end
    tick = medianFilter(tick, 3);
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

% Tinh MFCC cua 1 tin hieu
% Input:
%   y: mang gia tri cua tin hieu
%   Fs: tan so lay mau (Hz)
%   localSpeech: mang thoi gian thuc cua nguyen am (s)
%   N: so dac trung MFCC
% Output:
%   result: mang gia tri MFCC cua tin hieu
function result = mfccOf1File(y, Fs, localSpeech, N)
    time = (localSpeech(2) - localSpeech(1)) / 3;
    sTime = localSpeech(1) + time;
    eTime = localSpeech(1) + 2 * time;
    result = zeros(N, 1);
    S = mfcc(y(floor(sTime * Fs) : floor(eTime * Fs)), Fs, ...
                    'LogEnergy','Ignore', ...
                    'NumCoeffs', N);
    [count, ~] = size(S);
    for k = 1 : count
        for index = 1 : N
            result(index) = result(index) + S(k, index);
        end
    end
    result = result ./ count;
end

function result = fftOf1File(y, Fs, localSpeech, frameTime, frameShift, NFFT)
    result = zeros(1, NFFT / 2);
    time = (localSpeech(2) - localSpeech(1)) / 3;
    sTime = localSpeech(1) + time;
    eTime = localSpeech(1) + 2 * time;
    [frames, frameCount] = frameDivide(y(floor(sTime * Fs) : floor(eTime * Fs)), Fs, frameTime, frameShift);
    h = hamming(frameTime * Fs);
    for i = 1 : frameCount
        frame = frames(i, :)';
        frame = h.*frame;
        p = abs(fft(frame, NFFT));
        p = p(1:length(p)/2);
        p(2 : end - 1) = 2 * p(2 : end - 1);
        result(1, :) = result(1, :) + p';
    end
    result = result ./ frameCount;
end

% Tim vector dac trung de xac dinh nguyen am
% Input:
%   frameTime: do dai frame (s)
%   frameShift: do lech frame (s)
%   N: so dac trung MFCC
%   K: so cum kmean
%   path: duong dan chua file tin hieu
%   folder: ten folder chua file tin hieu
%   file: ten file tin hieu
% Ouput:
%   fingerprints: dau van tay dung de xac dinh nguyen am
function [fingerprints, fftDB] = getFingerprints(frameTime, frameShift, N, K, path, folder, file, NFFT)
    MFCC = zeros(length(folder), N);
    FFT = zeros(length(folder), NFFT / 2);
    fingerprints = zeros(length(file) * K, N);
    fftDB = zeros(length(file) * K, NFFT / 2);
    index = 1;
    for j = 1 : length(file) 
        for i = 1 : length(folder)
            filename = char(strcat(path, folder(i), '\',  file(j), '.wav'));
            [y, Fs] = audioread(filename);
            [localSpeech, ~] = voicedDetection(y, Fs, frameTime, frameShift);
            MFCC(i, :) = mfccOf1File(y, Fs, localSpeech, N)';
            FFT(i, :) = fftOf1File(y, Fs, localSpeech, frameTime, frameShift, NFFT);
        end
        rng(1);
        [~, C] = kmeans(MFCC, K);
        [~, D] = kmeans(FFT, K);
        fingerprints(index : index + K - 1, :) = C;
        fftDB(index : index + K - 1, :) = D;
        index = index + K;
    end
    filename = char(strcat(path, 'fingerprints.xlsx'));
    delete(filename);
    writetable(array2table(fingerprints), filename, "WriteVariableNames", false);
    filename = char(strcat(path, 'fftDB.xlsx'));
    delete(filename);
    writetable(array2table(fftDB), filename, "WriteVariableNames", false);
end

% Ve vector dac trung
% Input:
%   fingerprints: dau van tay dung de xac dinh nguyen am
%   K: so cum kmean
%   file: ten file tin hieu
function plotFigure(fingerprints, fftDB, K, file)
    [j, ~] = size(fingerprints);
    figure('Name', 'Fingerprints of 5 vowels');
    index = 1;
    for i = 1 : K : j
        subplot(5, 1, index);
        plot(fingerprints(i : i + K - 1, :)');
        title(char(file(index)));
        index = index + 1;
    end
    [j, ~] = size(fftDB);
    figure('Name', 'FFT of 5 vowels');
    index = 1;
    for i = 1 : K : j
        subplot(5, 1, index);
        plot(fftDB(i : i + K - 1, :)');
        title(char(file(index)));
        index = index + 1;
    end
end

% Tinh khoang cach eiclidean
% Input:
%   x: vector x
% y: vector y
% Output:
%   d: khoang cach cua 2 vector x va y
function d = euclidean(x, y)
    d = sum((x - y).^2).^0.5;
end

% Tinh do chinh xac cua thuat toan
% Input:
%   result: ket qua cua thuat toan
%   file: ten file tin hieu
%   folderKT: ten folder chua file tin hieu kiem thu 
% Output:
%   accuracy: do chinh xac cua thuat toan
function accuracy = calAccuracy(result, file, folderKT)
    count = 0;
    for i = 1 : length(result)
        for j = 1 : length(file)
            if (strcmp(result(i, j), file(j)) == true)
                count = count + 1;
            end
        end
    end
    accuracy = round(count / (length(file) * length(folderKT))  * 100, 2);
end
 
% Chay file tin hieu kiem thu
% Input:
%   fingerprints: dau van tay dung de xac dinh nguyen am
%   pathKT: duong dan chua file tin hieu kiem thu
%   folderKT: ten folder chua file tin hieu kiem thu 
%   file: ten file tin hieu
%   frameTime: do dai frame (s)
%   frameShift: do lech frame (s)
%   N: so dac trung MFCC
%   K: so cum kmean
% Output:
%   result: ket qua cua thuat toan
%   accuracy: do chinh xac cua thuat toan
function [result, accuracy, resultFFT, accuracyFFT] = runTesting(fingerprints, fftDB, pathKT, folderKT, file, frameTime, frameShift, N, K, NFFT)
    [count, ~] = size(fingerprints);
    t1 = zeros(1, count);
    t2 = zeros(1, count);
    result = cell(length(folderKT), length(file));
    resultFFT = cell(length(folderKT), length(file));
    for i = 1 : length(folderKT)
        for j = 1 : length(file)
            filename = char(strcat(pathKT, folderKT(i), '\',  file(j), '.wav'));
            [y, Fs] = audioread(filename);
            [localSpeech, ~] = voicedDetection(y, Fs, frameTime, frameShift);
            MFCC = mfccOf1File(y, Fs, localSpeech, N);
            FFT = fftOf1File(y, Fs, localSpeech, frameTime, frameShift, NFFT);
            for k = 1 : count
                t1(k) = euclidean(MFCC', fingerprints(k, :));
                t2(k) = euclidean(FFT, fftDB(k, :));
            end
            [~, index] = min(t1);
            result(i, j) = file(floor((index - 1) / K) + 1);
            [~, index] = min(t2);
            resultFFT(i, j) = file(floor((index - 1) / K) + 1);
        end
    end
     accuracy = calAccuracy(result, file, folderKT)
     accuracyFFT = calAccuracy(resultFFT, file, folderKT)
end







