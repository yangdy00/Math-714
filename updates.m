function [unew, ucurr] = updates(ucurr,uold, dt)
N = size(ucurr,1) - 1;
uL = Laplacian(ucurr);
uLL = Laplacian(uL);   %take Laplacian twice

unew = -uold + 2*ucurr + ((dt)^2)*uL + (1/12)*((dt)^4)*uLL;

%boundary conditions
unew(1,:) = 0;
unew(:,1) = 0;
unew(N+1,:) = 0;
unew(:,N+1) = 0;
end

