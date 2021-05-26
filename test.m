% Standard diffusion code, that solves
% du/dt = D*(d^2u/dx^2)
% using the forward Euler method.

D = 1.0; % Diffusion coefficient
Nx = 300; % Number of grid points in our simulation model
Dx = 0.1; % Spacing between grid points in our model, larger # makes larger system length

Dt = 0.1*(Dx*Dx)/D; % Timestep size (choose to be numerically stable)
Nt = 100000; % Number of timesteps to run

itplot = 100; % Plot every itplot timesteps

x = (0:(Nx-1))*Dx; % Define the coordinates of the gridpoints on the spatial grid

% % Initial conditions (u(x,t=0)) :
% % to simulate a forced firing, we will later make automatic
% % Start off the simulation with the dynamical variable u(x,t) equal
% % to a bell-shaped distribution in space, with center mu0, and
% % standard deviation sigma0:
% % mu0 = 0.5*(Nx-1)*Dx; % Put the center of the bell-shaped curve at the center of the system:
% mu0 = 0; % curve starts at left hand side
% % mu0 = 10; % curve starts at right hand side
% sigma0 = 1.0; % changes wideness and height of initial bump of wave
% u = exp(-(x-mu0).^2/(2*sigma0^2));
u = 0.01*rand(1,Nx);

% Define an auxiliary array, used in the calculation,
% which is the same size as the u array:
u_new = zeros(1,Nx);
v_new = zeros(1,Nx);

% PARAMETERS and ARRAY ALLOCATION
threshold = 0.1; % cell fires when u > threshold
magnify = 50; % makes firing term large enough
a = 0.8;
b_val = 0.05;
v = zeros(1, Nx); 
b = zeros(1,Nx); % define b for each cell 
for ix = 1:Nx % for each cell...
    % if we change this value to be smaller, we may not have enough mass to
    % generate a pulse and push it out of sinus node region
    if (ix<50) % if the cell is in the left region, composed of 49 cells
        b(ix) = -0.1; % set the value of b in this cell to avoid neg number
    else % right region
        b(ix) = b_val; % otherwise, set b i this cell to previous value 0.05
    end
end
N_branches = 3; % number of branches in fan-out used to model atrium
Nx_branch = 300; % number of cells in each branch
i_branch = 1; % which branch to plot
b_branch = 0.05; %b_value for branches
% extend coordinate system
x_branch = (Nx:(Nx + Nx_branch - 1)) * Dx;
% define arrays to keep track of u and v in each cell of each branch
u_branch = 0.01 * rand(N_branches, Nx_branch);
v_branch = zeros(N_branches, Nx_branch);
% these will hold new u and v values
u_branch_new = zeros(N_branches, Nx_branch);
v_branch_new = zeros(N_branches, Nx_branch);

% Timestep loop:
for it = 1:Nt
    
    for ix = 2:(Nx-1) % for all the interior points of the grid...
        % Advance the value of u to the next timestep, u_new, using
        % the forward Euler method in time to represent du/dt,
        % and centered differencing in space, to represent d^2/dx^2:
        term_1 = D * (u(ix-1) - 2*u(ix) + u(ix+1)) / Dx^2;
        % term_excite = magnify * u(ix) * (1 - u(ix)) * (u(ix) - threshold);
        threshold_2 = (v(ix) + b(ix)) / a;
        term_excite = magnify * u(ix) * (1 - u(ix)) * (u(ix) - threshold_2);
        % this models du/dt = D*(d^2u/dx^2) + 50*u*(1-u)*(u-0.1)
        % u_new(ix) = u(ix) + Dt * (term_1 + term_excite);
        % this models du/dt = D*(d^2u/dx^2) + 50*u*(1-u)*(u-0.1)
        temr_2 = 0;
        for i_branch = 1:N_branches
            term_2 = term_2 + D * (u(ix)-u_branch(i_branch,1)) / Dx^2;
        end
        u_new(ix) = u(ix) + Dt * (term_2 + term_excite);
    end
    for ix = 1:Nx
        v_new(ix) = v(ix) + Dt*(u(ix)-v(ix));
    end
    
    % update each branch
%     for i = 1:Nx
%         if i == 1
%             
%         else if i == Nx
%             for i_branch = 1:N_branches
%                 u_branch_new(i_branch, Nx_branch) = ...
%                     u_branch_new(i_branch, Nx_branch - 1);
%             end
%         else
%             
%         end
%     end
    
    u_branch = u_branch_new;
    v_branch = v_branch_new;
    
    % Enforce Neumann boundary conditions (du/dx=0) on the ends
    % of the system:
    u_new(1) = u_new(2);
    u_new(Nx) = u_new(Nx-1);
    
    % Set u = u_new for use in the next timestep:
    u = u_new;
    v = v_new;
    
    % for plotting
    x_plot = [x, x_branch];
    u_plot = [u, u_branch(i_branch_plot,:)];
    v_plot = [v, v_branch(i_branch_plot,:)];
    
    % Plot every so often:
    if (mod(it,itplot)==0)
        figure(1);
%         plot(x,u,x,v);
%         axis([x(1),x(Nx),0,1]); % define the plot axes
%         title(sprintf('u vs. x at time %f',it*Dt));
%         xlabel('x'); ylabel('u');
        plot(x_plot,u_plot,x_plot,v_plot);
        axis([x_plot(1),x_plot(end),0,1]); % define the plot axes
        title(sprintf('u and v vs. x at time %f',it*Dt));
        xlabel('x'); ylabel('u and v');
        drawnow;
    end
    
end