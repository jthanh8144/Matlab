n = -4:4;
x_n = 2*n;
E = sum(abs(x_n).^2);
display(E);
P = E/length(n);
display(P);