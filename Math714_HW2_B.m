error = 1;
N = 10;
while error > 10^-2
    N = N+1;
    x = 0:1/N:1;
    f = exp(-400*(x-0.5).^2);
    xq = 0:1/(1000*N):1;
    error = max(abs(interp1(x,f,xq,'linear')-exp(-400*(xq-0.5).^2)));
end
disp(N)
