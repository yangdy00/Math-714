% Script to benchmark the convergence of the numerical method

% number of refinements
Np = 4;
% arrays to store the error and mesh size
err = zeros(Np,1);
h = zeros(Np,1);

for idxRef = 1:Np

    p = 2^idxRef;
    n = 10*p; 
    m = 10*p;

    dx = 1/(n+1);
    dy = 1/(m+1);

    x = 0:dx:1;
    y = 0:dy:1;

    [X Y] = meshgrid(x,y) ;
    
    tol = 1e-6;
    [u,iter] = Jacobi_iterations(tol,n,m,dx,dy);

    % sampling the analytical solution 
    sol = cos(2*pi.*Y(:,2:end-1)).*(exp(-2*pi.*X(:,2:end-1))*exp(4*pi)...
        - exp(2*pi.*X(:,2:end-1)))/( exp(4*pi)-1 );
    u = u(2:n+1,:)';
    % computing the relative l^2 error
    err(idxRef) = max(abs(sol(:) - u(:)));
    % saving the step size
    h(idxRef) = dx;

end

% Plotting the results 
figure(1); clf();
loglog(h,err,'o-', 'LineWidth', 2)
hold on; 
loglog(h, h.^2, 'LineStyle', '-')

ax = gca;
ax.YAxis.FontSize = 13;
ax.XAxis.FontSize = 13;

title('Error', 'FontSize', 24);
xlabel('$h$','Interpreter','latex', 'FontSize', 24)
ylabel('relative $\ell^\infty$ error','Interpreter','latex', 'FontSize', 24)


lgd = legend("error", "$\mathcal{O}(h^2)$",'FontSize', 24,...
       'Interpreter','latex');
lgd.Location = 'northwest';