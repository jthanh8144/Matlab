clear;
frame_time = 0.02;
F0_Min = 70;
F0_Max = 450;
threshold = 0.3;
path = 'D:\Down\K� 5\XLTH\BT nhom\TinHieuHuanLuyen\';
file = {'phone_F1', 'phone_M1', 'studio_F1', 'studio_M1'};
for i = 1:4
    figure(i);
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    F0 = find_F0(y, Fs, frame_time, F0_Min, F0_Max, threshold);
    F0_t = medianFilter(F0, 1);
    
    subplot(4,1,3);
    plot (F0_t, '.');
    subplot(4,1,4);
    plot(y(1:end));
    [F0_mean, F0_std] = mean(F0_t)
    %RWfile(F0_mean, F0_std, char(file(i)));
    [delta_mean, delta_std] = delta(F0_mean, F0_std, path, file(i))
end

function yy = ACFunc(y, frame_length)
    yy = zeros(1, round(frame_length));
    for i = 0:length(y)-1
        yy(i+1)= 0;
        for j = 1:length(y)-i
            yy(i+1) = yy(i+1) + y(j)*y(i+j);
        end
    end
end

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

function F0 = find_F0(y, Fs, frame_time, F0_Min, F0_Max, threshold)
    N = length(y);
    frame_length = frame_time * Fs;
    frame_count = N / frame_length;
    F0 = zeros(floor(frame_count), 1);
    c1 = 0; c2 = 0;
    for i = 0:floor(frame_count)-1
        yt = y(i*frame_length+1:(i+1)*frame_length);
        yy = ACFunc(yt, frame_length);
        if(max(yt) > 0.01)
            F0(i+1) = find_F0_1frame(yy, Fs, F0_Min, F0_Max, threshold);
        end
        if (F0(i+1) == 0)
            if (c2 == 0) 
                subplot(4,1,2);
                plot (yy);
                c2=1;
            end
        else
            if (c1 == 0) 
                subplot(4,1,1);
                plot (yy);
                c1=1;
            end
        end
    end
end

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

function RWfile(F0_mean, F0_std, filename)
    path = 'D:\Down\K� 5\XLTH\BT nhom\test\';
    filename = char(strcat(path, filename, '.lab'));
    f = fopen(filename, 'a');
    fprintf(f, 'F0mean %.2f\n', F0_mean);
    fprintf(f, 'F0std %.2f\n', F0_std);
    fclose(f);
end

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
