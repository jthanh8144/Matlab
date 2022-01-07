clear;
frameTime = 0.02;
frameShift = 0.01;
path = 'D:\Down\Kì 5\XLTH\CK\TinHieuHuanLuyen\';
file = { '01MDA', '02FVA', '03MAB', '06FTB' };

for i = 1 : 4
    filename = char(strcat(path, file(i), '.wav'));
    [y, Fs] = audioread(filename);
    data = readFile(path, file(i));
    result = cal(data, y, Fs, frameTime, frameShift)
end

function data = frameDivide(y, Fs, frameTime, frameShift)
    N = length(y);
    frameLength = frameTime * Fs;
    shiftLength = frameShift * Fs;
    frameCount = floor((N - frameShift * Fs) / (frameLength - shiftLength));
    temp = 0;
    for i = 1:frameCount
        data(i, :) = y(temp + 1 : temp + frameLength);
        temp = temp + shiftLength;
    end
end

function [ste] = shortTimeEnergy(data)
    [frameCount, ~] = size(data);
    for i = 1:frameCount
        ste(i) = sum(data(i,:).^2);
    end
end

function data = readFile(path, name)
    filename = char(strcat(path, name, '.lab'));
    f = fopen(filename, 'r');
    line = {};
    i = 0;
    index = 0;
    while ~feof(f)
       i = i + 1;
       line(i) = cellstr(fgets(f));
    end
    for k = 1 : i - 2
        index = index + 1;
        temp= split(line(k));
        data(index) = str2double(temp(1));
        index = index + 1;
        data(index) = str2double(temp(2));
    end
end

function result = cal(data, y, Fs, frameTime, frameShift)
    index = 0;
    maxSTE = max(shortTimeEnergy(frameDivide(y, Fs, frameTime, frameShift)));
    
    for i = 1 : length(data) / 2
        index = index + 1;
        startPos = floor(data(index) * Fs + 1);
        index = index + 1;
        endPos = floor(data(index) * Fs);
        yy = y(startPos : endPos);
        yy_d = frameDivide(yy, Fs, frameTime, frameShift);
        ste = shortTimeEnergy(yy_d);
        result(i) = (sum(ste) / length(ste)) / maxSTE;
    end
end





