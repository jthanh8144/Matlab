clear;
path = 'D:\Down\Kì 5\XLTH\BT nhom\TinHieuHuanLuyen\';
file = {'phone_F1', 'phone_M1', 'studio_F1', 'studio_M1'};
i = 1;
frame_time = 0.02;
filename = char(strcat(path, file(i), '.wav'));
[y, Fs] = audioread(filename);
frame_length = frame_time * Fs;
a = 0.53;
b = 1.14;
yy = ACFunc(y(a*Fs:b*Fs), frame_length);

function yy = ACFunc(y, frame_length)
    yy = zeros(1, round(frame_length));
    for i = 0:length(y)-1
        yy(i+1)= 0;
        for j = 1:length(y)-i
            yy(i+1) = yy(i+1) + y(j)*y(i+j);
        end
    end
end