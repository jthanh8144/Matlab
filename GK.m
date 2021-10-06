% hzMin = 70 hzMax = 400
% chon nguong 30% cua gia tri cuc bo toan cuc r(0), thu voi 30 50 70
% [x, Fs] = audioread('D:\Down\Kì 5\XLTH\TinHieuMau\lab_female.wav');
%f, t = PitchTimeAmdf(x, 240, , Fs);

% khung; hzMin; hzMax; nguong, ChieuDaiCuaSoBoLoc là các thông s? m?c ??nh nghe
% khung = 0.02; hzMin = 75;
% hzMax = 350;
% nguong = 0.3;
% chieuDaiCuaSoBoLoc = 3;

% lab_female.wav, studio_female.wav --> hzMin = 120; hzMax = 400;

% disp('(1)lab_female (2)lab_male (3)studio_female (4)studio_male (other)finish');
% prompt = 'Chon file tin hieu:';
% i = input(prompt);


clear;
path = 'D:\Down\Kì 5\XLTH\TinHieuMau\';
file = ['lab_female.wav', 'lab_male.wav', 'phone_female.wav', 'phone_male.wav', 'studio_female.wav', 'studio_male.wav'];

% path + file(1)
[y,Fs] = audioread('D:\Down\Kì 5\XLTH\TinHieuMau\lab_male.wav');
khung = 0.02;
hzMin = 75;
hzMax = 350;
nguong = 0.3;
chieuDaiCuaSoBoLoc = 3;
%hzMin = 120; hzMax = 400;
[f0, x, yy] = AMDF(y, Fs, khung, hzMin, hzMax, nguong);
figure(2);
subplot(2,1,1);
stem(f0);
subplot(2,1,2);
plot(y);


