% chuong trinh lam tron 1 tin hieu de khu nhieu 
% dung bo loc lay trung binh cong 3 diem (3-points moving-averaging filter)
% co PTSP: y[n] = 1/3(x[n-1]+x[n]+x[n+1])

clf;                            % clear figures
L = 51;                         % do dai tin hieu
n = 0:L-1;                      % bien thoi gian roi rac
d = 1.5*randn(1,L);             % sinh tin hieu Gausian noise d[n] (0.5 la bien do nhieu)
s = 2*n.*(0.9.^n);              % sinh tin hieu goc s[n] = 2n(0.9)^n
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

% cach 1: dich thoi gian, lam tron tin hieu theo CT y[n] = 1/3(x[n-1]+x[n]+x[n+1])
x1 = [x(2:L), 0]        % x1[n] = x[n+1]
x2 = [x]                  % x2[n] = x[n]
x3 = [0, x(1:L-1)]    % x3[n] = x[n-1]

subplot(4,1,2)
% ve do thi x[n-1],x[n],x[n+1]
plot(n,x1,'r-.',n,x2,'b-.',n,x3,'k-.');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('x[n+1]','x[n]','x[n-1]');
title('time-shifted signals of x[n]');

y1 = 1/3*(x1+x2+x3)     % lay trung binh voi M=1
subplot(4,1,3)
% ve do thi y1[n] vs. s[n]
plot(n,y1(1:L),'r-',n,s(1:L),'b-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]','s[n]');
title('3-points smoothed y1[n] vs. original signal s[n]');

% cach 2: dung ham tinh tong chap conv()
% ghep noi tiep he som 1 don vi va he lay TB cong nhan qua
h = 1/3 * ones(1,3);    % h[n] = [1/3, 1/3, 1/3]
y2 = conv(x1, h);         % y2[n] = x1[n] * h[n]

% add
% cach 3: dung for
x0 = [0 x 0];
y3 = zeros(1, 51);
for i = 2 : (length(x0) - 1)
    y3(i-1) = (x0(i) + x0(i-1) + x0(i+1))/3;
end

% ve do thi y2[n] vs. s[n]
%hold on
subplot(4,1,4)
plot(n,y1,'r-',n,y2(1:L),'b-' , n, y3(1:L), 'g-');
xlabel('Chi so thoi gian n');
ylabel('Bien do');
legend('y1[n]','y2[n]','y3[n]');
title('cach1 vs. cach2 vs. cach3');



