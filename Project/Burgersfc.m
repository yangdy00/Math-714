function z = Burgersfc(u, range, ep, occupy)
% u_t+u*u_x=epsilon*u_xx
[v, w] = Derivative2fc(u, range, occupy);
z = ep * w - u.*v;

 
