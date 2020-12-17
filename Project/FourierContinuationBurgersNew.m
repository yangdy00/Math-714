dt = 0.01;
T = 25;
nt = round(T/dt);
ts = [0:dt:T]';
dx = 0.01;
xleft = -7;
xright = 7;
nx = round((xright-xleft)/dx);
xs = [xleft:dx:xright]';
ep = 1e-5;
u = zeros(nx+1,nt+1);
u(:,1) = 0.008*exp(-(0.2*xs).^10).*sin(2*pi*xs); % boundary condition
ufc = u;

% u_t+u_x=0
% Regular Fourier-based spectral method
u(:, 2) = u(:, 1) + Burgers(u(:, 1), xright-xleft, ep) * dt;
for i = 2:nt
    if mod(i, 1000) == 0
        disp(i)
    end
    u(:, i+1) = u(:, i-1) + 2 * Burgers(u(:, i), xright-xleft, ep) * dt;
end

% Fourier continuation method
range = 1;
nxrange = range/dx;
occupy = 0.8; % number of nodes overlapping at each end

step = Burgersfc(ufc(1:0.5*nxrange+6, 1), 0.5*range+6*dx, ep, occupy); % add 5 data points on each end
ufc(1:0.5*nxrange+1,2) = ufc(1:0.5*nxrange+1, 1) + step(1:end-5) * dt; % first section has 1 extra node at 0
for j = 1:(xright-xleft)/range-1
    step = Burgersfc(ufc((j-0.5)*nxrange-3:(j+0.5)*nxrange+6, 1), range+10*dx, ep, occupy);
    ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, 2) = ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, 1) + step(6:end-5) * dt;
end
step = Burgersfc(ufc((j+0.5)*nxrange-3:end, 1), 0.5*range+5*dx, ep, occupy);
ufc((j+0.5)*nxrange+2:end, 2) = ufc((j+0.5)*nxrange+2:end, 1) + step(6:end) * dt; % cannot extend to the right

for i = 2:nt
    if mod(i, 1000) == 0
        disp(i)
    end
    step = Burgersfc(ufc(1:0.5*nxrange+6, i), 0.5*range+6*dx, ep, occupy);
    ufc(1:0.5*nxrange+1,i+1) = ufc(1:0.5*nxrange+1, i-1) + 2 * step(1:end-5) * dt;
    for j = 1:(xright-xleft)/range-1
        step = Burgersfc(ufc((j-0.5)*nxrange-3:(j+0.5)*nxrange+6, i), range+10*dx, ep, occupy);
        ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, i+1) = ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, i-1) + 2 * step(6:end-5) * dt;
    end
    step = Burgersfc(ufc((j+0.5)*nxrange-3:end, i), 0.5*range+5*dx, ep, occupy);
    ufc((j+0.5)*nxrange+2:end, i+1) = ufc((j+0.5)*nxrange+2:end, i-1) + 2 * step(6:end) * dt;
    
%     points = linspace(ufc(0.5*nxrange+1, i+1), ufc(0.5*nxrange+7, i+1), 7);
%     ufc(0.5*nxrange+2:0.5*nxrange+6, i+1) = points(2:end-1);
%     for j = 1:(xright-xleft)/range-2
%         points = linspace(ufc((j+0.5)*nxrange-4, i+1), ufc((j+0.5)*nxrange+7, i+1), 12);
%         ufc((j+0.5)*nxrange-3:(j+0.5)*nxrange+6, i+1) = points(2:end-1);
%     end
%     points = linspace(ufc((j+1.5)*nxrange-4, i+1), ufc((j+1.5)*nxrange+2, i+1), 7);
%     ufc((j+1.5)*nxrange-3:(j+1.5)*nxrange+1, i+1) = points(2:end-1);
end

figure
subplot(3,2,1)
hold on
plot(xs,u(:,1))
plot(xs,ufc(:,1))
hold off
subplot(3,2,2)
hold on
plot(xs,u(:,5/dt+1))
plot(xs,ufc(:,5/dt+1))
hold off
subplot(3,2,3)
hold on
plot(xs,u(:,10/dt+1))
plot(xs,ufc(:,10/dt+1))
hold off
subplot(3,2,4)
hold on
plot(xs,u(:,15/dt+1))
plot(xs,ufc(:,15/dt+1))
hold off
subplot(3,2,5)
hold on
plot(xs,u(:,20/dt+1))
plot(xs,ufc(:,20/dt+1))
hold off
subplot(3,2,6)
hold on
plot(xs,u(:,25/dt+1))
plot(xs,ufc(:,25/dt+1))
hold off
