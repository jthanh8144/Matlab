clear;
frame_time = 0.02;
F0_Min = 70;
F0_Max = 450;
threshold = 0.3;
%path = 'D:\Down\Kì 5\XLTH\BT nhom\TinHieuHuanLuyen\';
%file = {'phone_F1', 'phone_M1', 'studio_F1', 'studio_M1'};
path = 'D:\Down\Kì 5\XLTH\BT nhom\TinHieuKiemThu\';
file = {'phone_F2', 'phone_M2', 'studio_F2', 'studio_M2'};
for i = 1:4
    figure('Name', char(file(i)));
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    F0 = find_F0(y, Fs, frame_time, F0_Min, F0_Max, threshold);
    F0_t = medianFilter(F0, 1);
    
    subplot(4,1,3);
    t = frame_time:frame_time:frame_time*length(F0_t);
    plot (t, F0_t, '.');
    xlabel('Time(s)');
    ylabel('F0 value(Hz)');
    title('Fundamental frequency (F0)');
    n=1:length(y);
    t = n/Fs;
    subplot(4,1,4);
    plot(t, y);
    xlabel('Time(s)');
    ylabel('Signal');
    title(['Input signal: ', char(file(i))]);
    
    [F0_mean, F0_std] = mean(F0_t)
    %RWfile(F0_mean, F0_std, path, char(file(i)));
    [delta_mean, delta_std] = delta(F0_mean, F0_std, path, file(i))
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

% Tim f0 cua mot khung
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

% Tim F0 cua tin hieu
% Input:
%   y: mang gia tri cua tin hieu
%   Fs: tan so lay mau (Hz)
%   frame_time: do dai cua mot khung (s)
%   F0_Min: tan so gioi han nho nhat
%   F0_Max: tan so gioi han lon nhat
%   threshold: nguong
% Output:
%   F0: mang gia tri F0 (Hz)
function F0 = find_F0(y, Fs, frame_time, F0_Min, F0_Max, threshold)
    c1 = 0; c2 = 0;
    N = length(y);
    frame_length = frame_time * Fs;
    frame_count = N / frame_length;
    F0 = zeros(floor(frame_count), 1);
    for i = 0:floor(frame_count)-1
        yt = y(i*frame_length+1:(i+1)*frame_length);
        yy = ACFunc(yt, frame_length);
        if(max(yt) > 0.01)
            F0(i+1) = find_F0_1frame(yy, Fs, F0_Min, F0_Max, threshold);
        end
        if (F0(i+1) == 0)
            if (c2 < 2)
                n=1:length(yy);
                t = n/Fs;
                subplot(4,1,2);
                plot(t, yy);
                xlabel('Time(s)');
                ylabel('ACF value');
                title('ACF value of 1 frame unvoice');
                c2 = c2 + 1;
            end
        else
            if (c1 < 2) 
                n=1:length(yy);
                t = n/Fs;
                subplot(4,1,1);
                plot(t, yy);
                xlabel('Time(s)');
                ylabel('ACF value');
                title('ACF value of 1 frame voice');
                c1 = c1 + 1;
            end
        end
    end
end

% Loc trung vi
% Input:
%   F0: mang gia tri F0 (Hz)
%   N: so phan tu lay o hai dau gia tri can tinh trung vi
% Output:
%   F0_med: mang gia tri F0 (Hz) da loc
function F0_med = medianFilter(F0, N)
    for i = 1:N
        F0 = [F0(1); F0; F0(length(F0))];
    end
    for i = N+1 : length(F0)-N
        t = [];
        for j = -N:N
            t = [t F0(i+j)];
        end
        t = sort(t);
        F0(i) = t(round((2*N + 1) / 2));
    end
    F0_med = F0(N+1:end-N);
end

% Ham tinh gia tri F0 trung binh va do lech chuan cua F0
% Input:
%   F0: mang gia tri F0 (Hz)
% Output:
%   F0_mean: gia tri F0 trung binh (Hz)
%   F0_std: gia tri do lech chuan cua F0 (Hz)
function [F0_mean, F0_std] = mean(F0)
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
end

% Doc du lieu va ghi ket qua F0_mean va F0_std vao file
% Input:
%   F0_mean: gia tri F0 trung binh (Hz)
%   F0_std: gia tri do lech chuan cua F0 (Hz)
%   filename: ten file du lieu
function RWfile(F0_mean, F0_std, path, filename)
    filename = char(strcat(path, filename, '.lab'));
    f = fopen(filename, 'a');
    fprintf(f, 'F0mean %.2f\n', F0_mean);
    fprintf(f, 'F0std %.2f\n', F0_std);
    fclose(f);
end

% Tinh do lech ket qua cua ham tu viet va file .lab in ra console
% Input:
%   F0_mean: gia tri F0 trung binh (Hz)
%   F0_std: gia tri do lech chuan cua F0 (Hz)
%   path: duong dan den file .lab
%   name: ten file .lab
% Output:
%   delta_mean: gia tri chenh lech cua F0 trung binh (Hz)
%   delta_std: gia tri chenh lech cua do lech chuan F0 (Hz)
function [delta_mean, delta_std] = delta(F0_mean, F0_std, path, name)
    filename = char(strcat(path, name, '.lab'));
    f = fopen(filename, 'r');
    line = {};
    i = 0;
    while ~feof(f)
       i = i + 1;
       line(i) = cellstr(fgets(f));
    end
    F0_mean_ = split(line(i-1));
    F0_mean_ = str2double(F0_mean_(2));
    F0_std_ = split(line(i));
    F0_std_ = str2double(F0_std_(2));
    delta_mean = abs(F0_mean - F0_mean_);
    delta_std = abs(F0_std - F0_std_);
    fclose(f);
end
