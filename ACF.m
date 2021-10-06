[y, Fs] = audioread('D:\Down\Kì 5\XLTH\TinHieuMau\phone_female.wav');
frame_time = 0.02;
F0_Min = 75;
F0_Max = 400;
threshold = 0.3;
F0 = find_F0(y, Fs, frame_time, F0_Min, F0_Max, threshold);
figure(2);
subplot(2,1,1);
stem (F0);
subplot(2,1,2);
plot(y(1:frame_time*Fs*length(F0)));

function yy = ACFunc(y, frame_length)
    yy = zeros(1, frame_length);
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
    f0 = Fs / t;
end

function F0 = find_F0(y, Fs, frame_time, F0_Min, F0_Max, threshold)
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
    end
end

