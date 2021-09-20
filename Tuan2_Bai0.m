n = -10:10;
u = [zeros(1, 10) 1 ones(1, 10)];

up1 = [u(2:length(u)) 1];
um5 = [zeros(1, 5) u(1:length(u) - 5)];
x1 = up1 - um5;

u0m = flip(u);
u2m = [ones(1, 2) u0m(1:length(u0m) - 2)];
x2 = n .* u2m;
x = x1 .* x2;

figure(1);
subplot(4, 2, 1);
stem(n, u, 'fill'); xlabel('Sample index n'); ylabel('u[n]');
subplot(4, 2, 2);
stem(n, up1, 'fill'); xlabel('Sample index n'); ylabel('u[n+1]');
subplot(4, 2, 3);
stem(n, um5, 'fill'); xlabel('Sample index n'); ylabel('u[n-5]');
subplot(4, 2, 4);
stem(n, x1, 'fill'); xlabel('Sample index n'); ylabel('x1[n]');
subplot(4, 2, 5);
stem(n, u0m, 'fill'); xlabel('Sample index n'); ylabel('u[-n]');
subplot(4, 2, 6);
stem(n, u2m, 'fill'); xlabel('Sample index n'); ylabel('u[2-n]');
subplot(4, 2, 7);
stem(n, x2, 'fill'); xlabel('Sample index n'); ylabel('x2[n]');
subplot(4, 2, 8);
stem(n, x, 'fill'); xlabel('Sample index n'); ylabel('x[n]');
