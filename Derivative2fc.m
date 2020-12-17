function [v, w] = Derivative2fc(u, range, occupy) % occupy < 1
N = floor((length(u)-1)/2);
indices = -N:N;
xs = linspace((1-occupy)/2, (1+occupy)/2, length(u))';
M = exp(1i*2*pi*xs*indices);
u_hat = M \ u;
v_hat = (-1i*2*pi*indices)'.* u_hat;
w_hat = (4*pi^2*indices.^2)'.* u_hat;
v = real(M*v_hat)*occupy/range;
w = real(M*w_hat)*occupy/range;
end
