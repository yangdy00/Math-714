function [u,iter] = Jacobi_iterations(tol,n,m,dx,dy)

% Start with the zero matrix
u0 = zeros(n+2,m+2);

% constant used
dxy2 = (dx^2*dy^2)/( 2*(dx^2+dy^2) );

for j = 1:m+2
    u0(1,j) = cos(2*pi*(j-1)*dy);
end

u=u0;
error = tol+1;
iter = 0;
while(error > tol)
    unew = u0;
    for i = 2:n+1
        unew(i,1)= dxy2*((u(i-1,1)+u(i+1,1))/(dx^2)+(2*u(i,2))/(dy^2)); % used ghost node
        for j = 2:m+1
            unew(i,j) = dxy2*((u(i-1,j)+u(i+1,j))/(dx^2)+(u(i,j-1)+u(i,j+1))/(dy^2));
        end
        unew(i,m+2)= dxy2*((u(i-1,m+2)+u(i+1,m+2))/(dx^2)+(2*u(i,m+1))/(dy^2));
    end
    error = norm(unew-u,2);
    u = unew;
    iter = iter + 1; 
end
end