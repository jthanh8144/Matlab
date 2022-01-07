clear;
% close all;
frameTime = 0.02;
frameShift = 0.01;
F0_Min = 70;
F0_Max = 400;
threshold_F0 = 0.3;
% path = 'D:\Down\Kì 5\XLTH\CK\TinHieuHuanLuyen\';
% file = { '01MDA', '02FVA', '03MAB', '06FTB' };
path = 'D:\Down\Kì 5\XLTH\CK\TinHieuKiemThu\';
file = { '30FTN', '42FQT', '44MTT', '45MDV' };

i = 1;
for i = 1 : 4
    figure('Name', char(file(i)));
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    [localSpeech, indexSpeech, index] = voicedDetection(y, Fs, frameTime, frameShift);
    yy = y(floor(localSpeech(9) * Fs) : floor(localSpeech(10) * Fs));
    [~, count] = frameDivide(y, Fs, frameTime, frameShift);
    [yyy_F0] = mientanso(yy, Fs);
    %yyy_F0 = medianFilter(yyy_F0, 3);
    %f0 = findF0(y, indexSpeech, Fs, frameTime, frameShift, 32768, F0_Max, F0_Min);
    subplot(4, 1, 1);
    plot(yyy_F0, '.');
    frameLength = frameTime * Fs;
    N = 2048; % 32768
    %shiftLength = frameShift * Fs;
    [frames, countF] = frameDivide(y, Fs, frameTime, frameShift);
    f0 = zeros(countF, 1);
    h = hamming(frameLength);
    data1 = frames(indexSpeech(1) : indexSpeech(1 + 1), :);
    for i = 1 : 2 : length(indexSpeech)
        startIndex = indexSpeech(i) - 1;
        data = frames(indexSpeech(i) : indexSpeech(i + 1), :);
        %data = h.*data;
        [dataCount, ~] = size(data);
        for j = 1 : dataCount
            frame = data(j, :)';
            frame = h.*frame;
            p = abs(fft(frame, N));
            p = p(1:length(p)/2+1);
            p(2 : end - 1) = 2 * p(2 : end - 1);
            freq = linspace(0, Fs/2, length(p));
            %if j == 1
                %[aa, bb] = findpeaks(p, freq);
            %end
            subplot(4, 1, 4);
            plot(p);
            if (findpeaks(frame, 'NPeaks', 1, 'SortStr', 'descend') > 0.03)
                [~, yPeak] = findpeaks(p, freq, 'MinPeakHeight', 2);%, 'MinPeakHeight', 1
                if (length(yPeak) > 2) 
                    f1 = yPeak(3) - yPeak(2);
                    f2 = yPeak(2) - yPeak(1);
                    if (f1 < F0_Max && f1 > F0_Min && f2 < F0_Max && f2 > F0_Min)
                        f0(startIndex + j) = (f1 + f2) / 2;
                    end
                end
            end
        end
    end
    t = (frameTime - frameShift) : (frameTime - frameShift) : (frameTime - frameShift) * length(f0);
    subplot(4, 1, 3);
    f0 = medianFilter(f0, 3);
    plot(t, f0, '.');
    countF0 = 0;
    sum = 0;
    for j = 1:length(f0)
        if (f0(j) > 0)
            sum = sum + f0(j);
            countF0 = countF0 + 1;
        end
    end
    mean = sum / countF0
end


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

function [ste] = shortTimeEnergy(data, frameCount)
    for i = 1:frameCount
        ste(i) = sum(data(i,:).^2);
    end
end

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

function [localSpeech, indexSpeech, index] = voicedDetection(y, Fs, frameTime, frameShift)
    [data, frameCount] = frameDivide(y, Fs, frameTime, frameShift);
    ste = shortTimeEnergy(data, frameCount);
    ste=ste./max(ste);
    threshold = findThreshold(ste, 10)
    t = (frameTime - frameShift) : (frameTime - frameShift) : (frameTime - frameShift) * length(ste);
    subplot(4, 1, 1);
    plot(t, ste,'r');
    title('Short time energy');
    xlabel('Time (s)');
    t = [0 : 1/Fs : length(y)/Fs];
    t = t(1 : end-1);
    subplot(4, 1, 2);
    plot(t, y);
    hold on;
    index = [];
    n = 0;
    for i = 1 : length(ste)
        if (ste(i) >= threshold)
            n = n + 1;
            index(n) = i;
        end
    end
    % tim diem thoi gian giao nhau giua tieng noi va khoang lang 
    indexSpeech(1) = index(1) - 1; 
    n = 2;
    for i=2:length(index)
        if ((frameTime - frameShift) * index(i) - (frameTime - frameShift) * index(i-1) > 0.3) 
            indexSpeech(n) = index(i-1);
            indexSpeech(n+1) = index(i) - 1;
            n=n+2;
        end
    end
    indexSpeech(n) = index(i);
    localSpeech = (frameTime - frameShift) * indexSpeech;
    % ve do thi phan doan tieng noi va khoang lang 
    x = -1:1; 
    for i = 1:length(localSpeech)
        if (rem(i,2) == 0)
            plot(localSpeech(i) * ones(size(x)), x,'r'); 
        else 
            plot(localSpeech(i) * ones(size(x)), x,'g'); 
        end
    end
    hold off; 
    legend('Signal', 'Voice start position', 'Voice end position'); 
    title('Voiced/Unvoiced detection'); 
    xlabel('Time (s)');
end

function f0 = findF0(y, localSpeech, Fs, frameTime, frameShift, N, F0_Max, F0_Min)
    h = hamming(frameTime * Fs);
    [~, count] = frameDivide(y(1 : localSpeech(1)), Fs, frameTime, frameShift);
    f0 = zeros(count, 1);
    for i = 1 : 2 : length(localSpeech)
        startIndex = length(f0);
        yy = y (floor(localSpeech(i)) : floor(localSpeech(i + 1)));
        [frames, frameCount] = frameDivide(yy, Fs, frameTime, frameShift);
        f0 = [f0 zeros(frameCount, 1)];
        for j = 1 : frameCount
            frame = frames(j, :)';
            frame = h.*frame;
            p = abs(fft(frame, N));
            p = p(1:length(p)/2+1);
            p(2 : end - 1) = 2 * p(2 : end - 1);
            
        end
    end



    frameLength = frameTime * Fs;
    %shiftLength = frameShift * Fs;
    frames = frameDivide(y, Fs, frameTime, frameShift);
    f0 = zeros(length(frames), 1);
    h = hamming(frameLength);
    for i = 1 : 2 : length(indexSpeech)
        startIndex = indexSpeech(i) - 1;
        data = frames(indexSpeech(i) : indexSpeech(i + 1));
        data = h.*data;
        for j = 1 : length(data)
            p = abs(fft(data(j), N));
            p = p(1:length(p)/2+1);
            p(2 : end - 1) = 2 * p(2 : end - 1);
            freq = linspace(0, Fs/2, length(p));
            if (findpeaks(data(j), 'NPeaks', 1, 'SortStr', 'descend') > 0.03)
                [~, yPeak] = findpeaks(p, freq, 'MinPeakHeight', 2);
                f1 = yPeak(3) - yPeak(2);
                f2 = yPeak(2) - yPeak(1);
                if (f1 < F0_Max && f1 > F0_Min && f2 < F0_Max && f2 > F0_Min)
                    f0(startIndex + j) = (f1 + f2) / 2;
                end
            end
        end
    end
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

function [yyy_F0] = mientanso(sig,fs)
    % The function is used for finding fundamental frequency of a signal in spectrum
    % sig is the input signal
    % fs is sampling frequency
    % The output is yyy_F0, which is the array containing calculated F0 in each windown
    %N = 32768;
    N = 1024;
    %N-point fft
    frame_len = round(0.02*fs); %Length of frame (30ms)
    half = round(frame_len/2); %For overlapping frame
    h = hamming(frame_len); %Hamming window function
    i=1; %Index of element in yyy_F0
    yyy_F0 = zeros(floor(length(sig)/half -1), 1);
    for k = 1 : length(sig)/half -1 %Vong lap cac frame
        range = (k-1)*half + 1:(frame_len + (k-1)*half); %Index of each element in window
        frame = h.*sig(range);
        %Value of each element in window
        %use FFT function (in Matlab) to analyze the spectrum of the frame
        P2 = abs(fft(frame,N));
        %The two-sided spectrum P2
        P1 = P2(1:length(P2)/2+1);
        %The single-sided spectrum P1
        P1(2:end-1) = 2*P1(2:end-1);
        freq=linspace(0,fs/2,length(P1)); %Spectrum of signal
        if findpeaks(frame,'NPeaks',1,'SortStr','descend') > 0.03 %Remove the noise by limiting the amplitude of signal in frame
            [y_value,y_peak] = findpeaks(P1,freq, 'MinPeakHeight', 1); %Find peaks (with the certain minimum Height to remove the noise in special cases) in the spectrum of signal
            %'MinPeakHeight',2
            z_a = y_peak(3) - y_peak(2); %Find 2 F0 by using 3 first peaks in spectrum
            z_b = y_peak(2) - y_peak(1);
            if z_a<400 && z_a>70
                if z_b<400 && z_b>70
                    yyy_F0(k) = (z_a+z_b)/2; %Final result is the mean of 2 F0
                    i=i+1;
                    %Put the result into the yyy_F0 and increase the index.
                end
            end
        end
    end
end




