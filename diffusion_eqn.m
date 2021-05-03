% Standard diffusion code, that solves
%
% du/dt = D*(d^2u/dx^2)
%
% using the forward Euler method.

D = 1.0; % Diffusion coefficient
Nx = 300; % Number of grid points in our simulation model
Dx = 0.1; % Spacing between grid points in our model

Dt = 0.1*(Dx*Dx)/D; % Timestep size (choose to be numerically stable)
Nt = 100000; % Number of timesteps to run

itplot = 100; % Plot every itplot timesteps

x = (0:(Nx-1))*Dx; % Define the coordinates of the gridpoints
% on the spatial grid

threshold = 0.2; %threshold value in u_new(ix)
magnify = 50;
a = 0.8;
b = zeros(1,Nx);
for ix = 1:Nx % for each cell
    if (ix<50) % if the cell is in the left region, composed of 49 cells
        b(ix) = -0.3; %set the value of b in this cell to avoid neg number
    else
        b(ix) = 0.05; %otherwise, set b i this cell to previous value 0.05
    end
end

u = 0.01*rand(1,Nx);

% Define an auxiliary array, used in the calculation,
% which is the same size as the u array:
u_new = zeros(1,Nx);
v = zeros(1,Nx);

% Timestep loop:
for it = 1:Nt
    
    for ix = 2:(Nx-1) % for all the interior points of the grid...
        % Advance the value of u to the next timestep, u_new, using
        % the forward Euler method in time to represent du/dt,
        % and centered differencing in space, to represent d^2/dx^2:
        
        term1 = D * (u(ix-1)-2*u(ix)+u(ix+1)) / Dx^2;
        term2 = magnify * u(ix) * (1-u(ix)) * (u(ix)-(v(ix)+b(ix))/a);
        u_new(ix) = u(ix) + Dt * (term1 + term2);
    end
    for ix = 1:Nx
        v_new(ix) = v(ix) + Dt*(u(ix)-v(ix));
    end
    
    % Enforce Neumann boundary conditions (du/dx=0) on the ends
    % of the system:
    u_new(1) = u_new(2);
    u_new(Nx) = u_new(Nx-1);
    
    % Set u = u_new for use in the next timestep:
    u = u_new;
    v = v_new;
    
    % Plot every so often:
    if (mod(it,itplot)==0)
        figure(1);
        plot(x,u,x,v);
        axis([x(1),x(Nx),0,1]); % define the plot axes
        title(sprintf('u vs. x at time %f',it*Dt));
        xlabel('x'); ylabel('u');
        drawnow;
    end
    
end