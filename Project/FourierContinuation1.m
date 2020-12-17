dt = 0.001;
T = 30;
nt = round(T/dt);
ts = [0:dt:T]';
dx = 0.01;
xleft = -10;
xright = 40;
nx = round((xright-xleft)/dx);
xs = [xleft:dx:xright]';
u = zeros(nx+1,nt+1);
u(:,1) = exp(-(0.2*xs).^10).*sin(2*pi*xs); % boundary condition
ufc = u;

% u_t+u_x=0
% Regular Fourier-based spectral method
u(:, 2) = u(:, 1) - Derivativefft(u(:, 1), xright-xleft) * dt;
for i = 2:nt
    if mod(i, 1000) == 0
        disp(i)
    end
    u(:, i+1) = u(:, i-1) - 2 * Derivativefft(u(:, i), xright-xleft) * dt;
end

% Fourier continuation method
range = 1;
nxrange = range/dx;
fringe = nxrange * 0.5; % number of nodes overlapping at each end

% % split through the middle
% deriv = Derivativefft(ufc(1:nxrange+fringe+1, 1), range+(fringe+1)*dx);
% ufc(1:nxrange+1,2) = ufc(1:nxrange+1, 1) - deriv(1:nxrange+1) * dt; % first section has 1 extra node at 0
% for j = 2:(xright-xleft)/range-1
%     deriv = Derivativefft(ufc((j-1)*nxrange-fringe+2:j*nxrange+fringe+1, 1), range+2*fringe*dx);
%     ufc((j-1)*nxrange+2:j*nxrange+1, 2) = ufc((j-1)*nxrange+2:j*nxrange+1, 1) - deriv(fringe+1:fringe+nxrange) * dt;
% end
% deriv = Derivativefft(ufc(j*nxrange-fringe+2:end, 1), range+fringe*dx);
% ufc(j*nxrange+2:end, 2) = ufc(j*nxrange+2:end, 1) - deriv(fringe+1:end) * dt; % cannot extend to the right
% 
% for i = 2:nt
%     deriv = Derivativefft(ufc(1:nxrange+fringe+1, i), range+(fringe+1)*dx);
%     ufc(1:nxrange+1,i+1) = ufc(1:nxrange+1, i-1) - 2 * deriv(1:nxrange+1) * dt;
%     if mod(i, 1000) == 0
%         disp(i)
%     end
%     for j = 2:(xright-xleft)/range-1
%         deriv = Derivativefft(ufc((j-1)*nxrange-fringe+2:j*nxrange+fringe+1, i), range+2*fringe*dx);
%         ufc((j-1)*nxrange+2:j*nxrange+1, i+1) = ufc((j-1)*nxrange+2:j*nxrange+1, i-1) - 2 * deriv(fringe+1:fringe+nxrange) * dt;
%     end
%     deriv = Derivativefft(ufc(j*nxrange-fringe+2:end, i), range+fringe*dx);
%     ufc(j*nxrange+2:end, i+1) = ufc(j*nxrange+2:end, i-1) - 2 * deriv(fringe+1:end) * dt;
% end

% % preserve middle
% deriv = Derivativefft(ufc(1:0.5*nxrange+fringe+1, 1), 0.5*range+(fringe+1)*dx);
% ufc(1:0.5*nxrange+1,2) = ufc(1:0.5*nxrange+1, 1) - deriv(1:0.5*nxrange+1) * dt; % first section has 1 extra node at 0
% % deriv = Burgers(ufc(1:1.5*nxrange+fringe+1, 1), 1.5*range+(fringe+1)*dx, ep);
% % ufc(0.5*nxrange+2:1.5*nxrange+1, 2) = ufc(0.5*nxrange+2:1.5*nxrange+1, 1) + deriv(0.5*nxrange+2:1.5*nxrange+1) * dt;
% for j = 2:(xright-xleft)/range-2
%     deriv = Derivativefft(ufc((j-0.5)*nxrange-fringe+2:(j+0.5)*nxrange+fringe+1, 1), range+2*fringe*dx);
%     ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, 2) = ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, 1) - deriv(fringe+1:fringe+nxrange) * dt;
% end
% % deriv = Burgers(ufc((j+0.5)*nxrange-fringe+2:end, 1), 1.5*range+fringe*dx, ep);
% % ufc((j+0.5)*nxrange+2:(j+1.5)*nxrange+1, 2) = ufc((j+0.5)*nxrange+2:(j+1.5)*nxrange+1, 1) + deriv(fringe+1:fringe+nxrange) * dt;
% deriv = Derivativefft(ufc((j+0.5)*nxrange-fringe+2:end, 1), 0.5*range+fringe*dx);
% ufc((j+0.5)*nxrange+2:end, 2) = ufc((j+0.5)*nxrange+2:end, 1) - deriv(fringe+1:end) * dt; % cannot extend to the right
% 
% for i = 2:nt
%     if mod(i, 1000) == 0
%         disp(i)
%     end
%     deriv = Derivativefft(ufc(1:0.5*nxrange+fringe+1, i), 0.5*range+(fringe+1)*dx);
%     ufc(1:0.5*nxrange+1,i+1) = ufc(1:0.5*nxrange+1, i-1) - 2 * deriv(1:0.5*nxrange+1) * dt;
% %     deriv = Burgers(ufc(1:1.5*nxrange+fringe+1, i), 1.5*range+(fringe+1)*dx, ep);
% %     ufc(0.5*nxrange+2:1.5*nxrange+1, i+1) = ufc(0.5*nxrange+2:1.5*nxrange+1, i-1) + 2 * deriv(0.5*nxrange+2:1.5*nxrange+1) * dt;
%     for j = 1:(xright-xleft)/range-1
%         deriv = Derivativefft(ufc((j-0.5)*nxrange-fringe+2:(j+0.5)*nxrange+fringe+1, i), range+2*fringe*dx);
%         ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, i+1) = ufc((j-0.5)*nxrange+2:(j+0.5)*nxrange+1, i-1) - 2 * deriv(fringe+1:fringe+nxrange) * dt;
%     end
% %     deriv = Burgers(ufc((j+0.5)*nxrange-fringe+2:end, i), 1.5*range+fringe*dx, ep);
% %     ufc((j+0.5)*nxrange+2:(j+1.5)*nxrange+1, i+1) = ufc((j+0.5)*nxrange+2:(j+1.5)*nxrange+1, i-1) + 2 * deriv(fringe+1:fringe+nxrange) * dt;
%     deriv = Derivativefft(ufc((j+0.5)*nxrange-fringe+2:end, i), 0.5*range+fringe*dx);
%     ufc((j+0.5)*nxrange+2:end, i+1) = ufc((j+0.5)*nxrange+2:end, i-1) - 2 * deriv(fringe+1:end) * dt;
% end

% % 1 step forward
% for i = 2:nt
%     deriv = Derivativefft(ufc(1:nxrange+fringe+1, i), range+(fringe+1)*dx);
%     ufc(1:nxrange+1,i+1) = ufc(1:nxrange+1, i) - deriv(1:nxrange+1) * dt;
%     if mod(i, 1000) == 0
%         disp(i)
%     end
%     for j = 2:(xright-xleft)/range-1
%         deriv = Derivativefft(ufc((j-1)*nxrange-fringe+2:j*nxrange+fringe+1, i), range+2*fringe*dx);
%         ufc((j-1)*nxrange+2:j*nxrange+1, i+1) = ufc((j-1)*nxrange+2:j*nxrange+1, i) - deriv(fringe+1:fringe+nxrange) * dt;
%     end
%     deriv = Derivativefft(ufc(j*nxrange-fringe+2:end, i), range+fringe*dx);
%     ufc(j*nxrange+2:end, i+1) = ufc(j*nxrange+2:end, i) - deriv(fringe+1:end) * dt;
% end

u_true = zeros(nx+1,nt+1);
for i = 1:T/dt+1
    u_true(:,i) = exp(-(0.2*(xs-ts(i))).^10).*sin(2*pi*(xs-ts(i)));
end

error = zeros(1, T/dt+1);
for i = 1:T/dt+1
    error(i) = norm(abs(u(:,i)-u_true(:,i)));
end

figure
subplot(2,2,1)
plot(xs,u(:,1));
subplot(2,2,2)
plot(xs,u(:,10/dt+1))
subplot(2,2,3)
plot(xs,u(:,20/dt+1))
subplot(2,2,4)
plot(xs,u(:,30/dt+1))

% figure
% subplot(2,2,1)
% plot(xs,ufc(:,1));
% subplot(2,2,2)
% plot(xs,ufc(:,10/dt+1))
% subplot(2,2,3)
% plot(xs,ufc(:,20/dt+1))
% subplot(2,2,4)
% plot(xs,ufc(:,30/dt+1))
