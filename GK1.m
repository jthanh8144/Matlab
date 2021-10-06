clear;
[y, fs] = audioread('D:\Down\Kì 5\XLTH\TinHieuMau\test\LA001.wav');
a = y(4000:5000);
sum = zeros(1,320);
for k = 320
    sum(k) = 0;
end
for k=1:320
    for i=1:32
        sum(k) = sum(k) + (a(i)*a(i+k));
        sum(k) = sum(k)/32;
    end
end
[pks, locs] = findpeaks(sum);
min(fs./diff(locs)), mean(fs./diff(locs)), max(fs./diff(locs));
disp(locs);
[mm, peak1_ind] = min((fs./diff(locs)));
period = locs(peak1_ind + 1) - locs(peak1_ind);
pitch_Hz = fs/period;
%disp(pitch_Hz);