A = 1;
F0 = 1000;
pha = 0;
t1 = 0;
t2 = 3;

Fs1 = 3*F0;
Ts1 = 1/Fs1;
n1 = t1:Ts1:t2;
x1 = A*cos(2*pi*F0*n1 + pha);

Fs2 = 1.5*F0;
Ts2 = 1/Fs2;
n2 = t1:Ts2:t2;
x2 = A*cos(2*pi*F0*n2 + pha);

figure(1);
subplot(2, 1, 1);
plot(n1(1:10), x1(1:10)); xlabel('Time'); ylabel('Audio Signal x1');
subplot(2, 1, 2);
plot(n2(1:10), x2(1:10)); xlabel('Time'); ylabel('Audio Signal x2');

%sound(x1, Fs1);
%sound(x2, Fs2);