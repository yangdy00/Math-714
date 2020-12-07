N = 100;
thetas = (1/N)*pi*(0:N);
xs = cos(thetas);
xs = (1/2)*xs + 0.5; % change it to between 0 and 1
ys = xs';
dt = 6*(1/N)^2; % smaller than the CFL condition

fineN = 600; % fine grid
finethetas = (1/fineN)*pi*(0:fineN);
finexs = cos(finethetas);
finexs = (1/2)*finexs + 0.5;
fineys = finexs';
finedt = 6*(1/fineN)^2;

uold = zeros(N+1,N+1);
ucurr = uold + Fx(ys)*Fx(xs)*dt + (1/2)*(dt^2)*Laplacian(uold) + (1/6)*(dt^3)*Laplacian(Fx(ys)*Fx(xs));
ucurr(1,:) = 0;
ucurr(N+1,:) = 0;
ucurr(:,1) = 0;
ucurr(:,N+1) = 0;

fineuold = zeros(fineN+1,fineN+1);
fineucurr = fineuold + Fx(fineys)*Fx(finexs)*finedt + (1/2)*(finedt^2)*Laplacian(fineuold) + (1/6)*(finedt^3)*Laplacian(Fx(fineys)*Fx(finexs));
fineucurr(1,:) = 0;
fineucurr(fineN+1,:) = 0;
fineucurr(:,1) = 0;
fineucurr(:,fineN+1) = 0;

fineRounds = round(dt / finedt);


for tmp=2:fineRounds
    [fineucurr, fineuold] = updates(fineucurr, fineuold, finedt);
end

Rounds = round(0.1/dt);
error = max(abs(fineucurr(1:diff:fineN+1,1:diff:fineN+1)-ucurr), [], 'all');

for rounds=2:Rounds

    if mod(rounds, 10) == 0
        disp(['Round ', num2str(rounds), ' out of ', num2str(Rounds)])
        disp(['Current error: ', num2str(error)]);
    end
    [ucurr, uold] = updates(ucurr, uold, dt);
    
    for tmp=1:fineRounds
       %if mod(tmp, 10) == 0
       %    disp(['findIteration ', num2str(tmp), ' out of ', num2str(fineRounds)]);
       %end

       [fineucurr, fineuold] = updates(fineucurr, fineuold, finedt);
    end

    diff = fineN/N;
    error = max(error, max(abs(fineucurr(1:diff:fineN+1,1:diff:fineN+1)-ucurr), [], 'all'));
    
end

maxE = zeros(1,3);
maxE(1) = error; % record

N = 200;
thetas = (1/N)*pi*(0:N);
xs = cos(thetas);
xs = (1/2)*xs + 0.5; % change it to between 0 and 1
ys = xs';
dt = 6*(1/N)^2; % smaller than the CFL condition

fineN = 600; % fine grid
finethetas = (1/fineN)*pi*(0:fineN);
finexs = cos(finethetas);
finexs = (1/2)*finexs + 0.5;
fineys = finexs';
finedt = 6*(1/fineN)^2;

uold = zeros(N+1,N+1);
ucurr = uold + Fx(ys)*Fx(xs)*dt + (1/2)*(dt^2)*Laplacian(uold) + (1/6)*(dt^3)*Laplacian(Fx(ys)*Fx(xs));
ucurr(1,:) = 0;
ucurr(N+1,:) = 0;
ucurr(:,1) = 0;
ucurr(:,N+1) = 0;

fineuold = zeros(fineN+1,fineN+1);
fineucurr = fineuold + Fx(fineys)*Fx(finexs)*finedt + (1/2)*(finedt^2)*Laplacian(fineuold) + (1/6)*(finedt^3)*Laplacian(Fx(fineys)*Fx(finexs));
fineucurr(1,:) = 0;
fineucurr(fineN+1,:) = 0;
fineucurr(:,1) = 0;
fineucurr(:,fineN+1) = 0;

fineRounds = round(dt / finedt);


for tmp=2:fineRounds
    [fineucurr, fineuold] = updates(fineucurr, fineuold, finedt);
end

Rounds = round(0.1/dt);
error = max(abs(fineucurr(1:diff:fineN+1,1:diff:fineN+1)-ucurr), [], 'all');

for rounds=2:Rounds

    if mod(rounds, 10) == 0
        disp(['Round ', num2str(rounds), ' out of ', num2str(Rounds)])
        disp(['Current error: ', num2str(error)]);
    end
    [ucurr, uold] = updates(ucurr, uold, dt);
    
    for tmp=1:fineRounds
       %if mod(tmp, 10) == 0
       %    disp(['findIteration ', num2str(tmp), ' out of ', num2str(fineRounds)]);
       %end

       [fineucurr, fineuold] = updates(fineucurr, fineuold, finedt);
    end

    diff = fineN/N;
    error = max(error, max(abs(fineucurr(1:diff:fineN+1,1:diff:fineN+1)-ucurr), [], 'all'));
    
end

maxE(2) = error; % record

N = 300;
thetas = (1/N)*pi*(0:N);
xs = cos(thetas);
xs = (1/2)*xs + 0.5; % change it to between 0 and 1
ys = xs';
dt = 6*(1/N)^2; % smaller than the CFL condition

fineN = 600; % fine grid
finethetas = (1/fineN)*pi*(0:fineN);
finexs = cos(finethetas);
finexs = (1/2)*finexs + 0.5;
fineys = finexs';
finedt = 6*(1/fineN)^2;

uold = zeros(N+1,N+1);
ucurr = uold + Fx(ys)*Fx(xs)*dt + (1/2)*(dt^2)*Laplacian(uold) + (1/6)*(dt^3)*Laplacian(Fx(ys)*Fx(xs));
ucurr(1,:) = 0;
ucurr(N+1,:) = 0;
ucurr(:,1) = 0;
ucurr(:,N+1) = 0;

fineuold = zeros(fineN+1,fineN+1);
fineucurr = fineuold + Fx(fineys)*Fx(finexs)*finedt + (1/2)*(finedt^2)*Laplacian(fineuold) + (1/6)*(finedt^3)*Laplacian(Fx(fineys)*Fx(finexs));
fineucurr(1,:) = 0;
fineucurr(fineN+1,:) = 0;
fineucurr(:,1) = 0;
fineucurr(:,fineN+1) = 0;

fineRounds = round(dt / finedt);


for tmp=2:fineRounds
    [fineucurr, fineuold] = updates(fineucurr, fineuold, finedt);
end

Rounds = round(0.1/dt);
error = max(abs(fineucurr(1:diff:fineN+1,1:diff:fineN+1)-ucurr), [], 'all');

for rounds=2:Rounds

    if mod(rounds, 10) == 0
        disp(['Round ', num2str(rounds), ' out of ', num2str(Rounds)])
        disp(['Current error: ', num2str(error)]);
    end
    [ucurr, uold] = updates(ucurr, uold, dt);
    
    for tmp=1:fineRounds
       %if mod(tmp, 10) == 0
       %    disp(['findIteration ', num2str(tmp), ' out of ', num2str(fineRounds)]);
       %end

       [fineucurr, fineuold] = updates(fineucurr, fineuold, finedt);
    end

    diff = fineN/N;
    error = max(error, max(abs(fineucurr(1:diff:fineN+1,1:diff:fineN+1)-ucurr), [], 'all'));
    
end

maxE(3) = error; % record



N= [ 100, 200, 300 ];
maxE = [ 0.0015925, 0.00075564, 0.00051251];


figure(1); clf();
loglog(N,maxE,'o-', 'LineWidth', 2, 'color', 'b')
hold on; 
loglog(N, (1./N).^4, 'LineStyle', '-', 'color', 'r')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error', 'FontSize', 24)
xlabel('N', 'FontSize', 18)
ylabel('max error relative to a fine grid', 'FontSize', 18)

function [y] = Fx(x)
y = exp(-400*((0.5*x+0.5)-0.5).^2);
end

