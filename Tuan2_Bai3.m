[x, Fs] = audioread('D:\Down\K� 5\XLTH\TinHieuMau\lab_female.wav');
sound_duration = length(x)/Fs;
t = 0:1/Fs:sound_duration;
t = t(1:end-1);
plot(t, x);
xlabel('Time');
ylabel('Audio Signal');
sound(x, Fs);
sound(x, Fs/2);
sound(x, Fs*2);