% tim x[2n] biet x[n] co dang xung tam giac
NL = 20;
NR = 10;
minA = 0;
maxA = 5;
n = -NL:NR;
x = [linspace(minA, maxA, NL) maxA linspace(maxA, minA, NR)];

x2 = zeros(1, length(x));
for i = 0:NL/2
    x2(NL + 1 - i) = x(NL + 1 - 2*i);
end
for i = 0:NR/2
    x2(NL + 1 + i) = x(NL + 1 + 2*i);
end

subplot(2, 1, 1);
stem(n, x, 'fill'); xlabel('n'); ylabel('x[n]');
subplot(2, 1, 2);
stem(n, x2, 'fill'); xlabel('n'); ylabel('x[2n]');