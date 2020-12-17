function w = Burgers(u, range, ep)
% u_t+u*u_x=epsilon*u_xx
ux = Derivativefft(u, range);
w = ep * Derivativefft(ux, range) - u.*ux;