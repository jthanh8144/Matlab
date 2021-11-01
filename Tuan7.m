% chuong trinh lam tron 1 tin hieu de khu nhieu 
% dung bo loc lay trung binh cong 3 diem (3-points moving-averaging filter)
% co PTSP: y[n] = 1/3(x[n]+x[n-1]+x[n-2])

clf;                            % clear figures
L = 51;                         % do dai tin hieu
n = 0:L-1;                      % bien thoi gian roi rac
d = 0.5*randn(1,L);             % sinh tin hieu Gausian noise d[n] (0.5 la bien do nhieu)
s = 5 * 0.9.^n;                 % sinh tin hieu goc s[n] = 2n(0.9)^n
x = s + d;                      % tin hieu co nhieu x[n]=s[n]+d[n]

% problem: co tin hieu bi nhieu x[n], di khoi phuc lai tin hieu goc s[n]
figure(1)                       % tao figure 1
%hold on
subplot(4,1,1)                  % tao luoi 3x1 plots va ve plot o vi tri #1
plot(n,d,'r-',n,s,'k--',n,x,'b-.'); % ve do thi d[n],s[n],x[n]
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('d[n]','s[n]','x[n]');
title('Noise d[n] vs. original s[n] vs. noisy signals x[n]');

% cach 1: dich thoi gian, lam tron tin hieu theo CT y[n] = 1/3(x[n]+x[n-1]+x[n-2])
x1 = [x];                         % x1[n] = x[n]
x2 = [0, x(1:L-1)];               % x2[n] = x[n-1]
x3 = [zeros(1, 2), x(1:L-2)];     % x3[n] = x[n-2]

subplot(4,1,2)
% ve do thi x[n],x[n-1],x[n-2]
plot(n,x1,'r-.',n,x2,'b-.',n,x3,'k-.');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('x[n]','x[n-1]','x[n-2]');
title('time-shifted signals of x[n]');

y1 = 1/3*(x1+x2+x3);     % lay trung binh voi M=1
subplot(4,1,3)
% ve do thi y1[n] vs. s[n]
plot(n,y1(1:L),'r-',n,s(1:L),'b-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]','s[n]');
title('3-points smoothed y1[n] vs. original signal s[n]');

% cach 2: dung for
x0 = [0 0 x];
y2 = zeros(1, 51);
for i = 3:length(x0)
    y2(i-2) = (x0(i) + x0(i-1) + x0(i-2))/3;
end

% cach 3: dung ham tinh tong chap conv()
% ghep noi tiep he som 1 don vi va he lay TB cong nhan qua
h = 1/3 * ones(1,3);    % h[n] = [1/3, 1/3, 1/3]
y3 = conv(x1, h);         % y2[n] = x1[n] * h[n]

% cach 4: dung ham cai dat bo loc filter()???
windowSize = 3;
b = (1 / windowSize) * ones(1, windowSize);
a = 1;
y4 = filter(b, a, x);

% ve do thi y2[n], y3[n] vs. s[n]
subplot(4,1,4)
plot(n, y1, 'r-', n, y2(1:L), 'b-', n, y3(1:L), 'm-', n, y4(1:L), 'g-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]', 'y2[n]', 'y3[n]', 'y4[n]');
title('cach1  vs. cach2 vs. cach3 vs. cach4');
