clear;
close all;
path = 'D:\Down\Kì 5\XLTH\ThucHanh\TH1\NguyenAmHuanLuyen-16k\';
folder = { '15MMH', '16FTH', '17MTH', '18MNK', '19MXK', '20MVK', '21MTL', '22MHL' };
file = { 'a', 'e', 'i', 'o', 'u' };

for i = 1 : length(folder)
    figure('Name', char(folder(i)));
    for j = 1 : length(file)
        filename = char(strcat(path, folder(i), '\',  file(j), '.wav'));
        [y, Fs] = audioread(filename);
        
        subplot(5, 1, j);
        spectrogram(y, floor(5*10^(-3)*Fs), floor(2*10^(-3)*Fs), 2048, Fs, 'yaxis');
        xlabel('Time(s)');
        title(strcat('Wideband vowel "', file(j), '"'));
    end
end

