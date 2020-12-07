function [L] = Laplacian(u)
N = size(u,1) - 1;
thetas = (1/N)*pi*(0:N);
x = cos(thetas);
y = x';

uxx = zeros(N+1,N+1); uyy = zeros(N+1,N+1);
for i = 1:N+1 % 2nd derivatives wrt x in each row
    v = u(i,:); 
    w = chebfft(v');
    w2 = chebfft(w);
    uxx(i,:) = w2';
end
for j = 1:N+1 % 2nd derivatives wrt y in each column
    v = u(:,j); 
    w = chebfft(v);
    w2 = chebfft(w);
    uyy(:,j) = w2;
end
L = uxx + uyy;
end

