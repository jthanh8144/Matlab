clear;
xn = audioread('D:\Down\K� 5\XLTH\TinHieuMau\test\LA001.wav');
E = sum(abs(xn).^2);
display(E);
P = E / length(xn);
display(P);