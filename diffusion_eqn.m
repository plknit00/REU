% Standard diffusion code, that solves
% du/dt = D*(d^2u/dx^2) using the forward Euler method.

% PLOT VALUES
D = 1.0; % Diffusion coefficient
Nx = 300; % Number of grid points in our simulation model
Dx = 0.1; % Spacing between grid points in our model, larger # makes larger system length
Dt = 0.1*(Dx*Dx)/D; % Timestep size (choose to be numerically stable)
Nt = 100000; % Number of timesteps to run
itplot = 100; % Plot every itplot timesteps
x = (0:(Nx-1))*Dx; % Define the coordinates of the gridpoints on the spatial grid
u = 0.01*rand(1,Nx);

% PARAMETERS and ARRAY ALLOCATION
u_new = zeros(1,Nx);
v_new = zeros(1,Nx);
magnify = 50; % makes firing term large enough
a = 0.8;
b_val = 0.05;
v = zeros(1, Nx); 
b = zeros(1,Nx); % define b for each cell 
for ix = 1:Nx 
    % if we change this value to be smaller, we may not have enough mass to
    % generate a pulse and push it out of sinus node region
    if (ix<50) % if the cell is in the left region, composed of 49 cells
        b(ix) = -0.1; % set the value of b in this cell to avoid neg number
    else % right region
        b(ix) = b_val; % otherwise, set b i this cell to previous value 0.05
    end
end

N_branches = 1; % number of branches in fan-out used to model atrium
Nx_branch = 300; % number of cells in each branch
i_branch_plot = 1; % which branch to plot
b_branch = 0.05; % b_value for branches
x_branch = (Nx:(Nx + Nx_branch - 1)) * Dx; % extend coordinate system
% define arrays to keep track of u and v in each cell of each branch
u_branch = 0.01 * rand(N_branches, Nx_branch);
v_branch = zeros(N_branches, Nx_branch);
u_branch_new = zeros(N_branches, Nx_branch);
v_branch_new = zeros(N_branches, Nx_branch);

% Timestep loop:
for it = 1:Nt
    for ix = 2:(Nx-1) % for all the interior points of the grid...
        % Advance the value of u to the next timestep, u_new, using
        % the forward Euler method in time to represent du/dt,
        % and centered differencing in space, to represent d^2/dx^2:
        threshold = (v(ix) + b(ix)) / a; % cell fires when u > threshold
        term_excite = magnify * u(ix) * (1 - u(ix)) * (u(ix) - threshold);
        term_2 = 0;
        for i_branch = 1:N_branches
            term_2 = term_2 + D * (u(ix)-u_branch(i_branch,1)) / Dx^2;
        end
        u_new(ix) = u(ix) + Dt * (term_2 + term_excite);
    end
    for ix = 1:Nx
        v_new(ix) = v(ix) + Dt*(u(ix)-v(ix));
    end
    
    % update each branch
    for i = 1:N_branches
        if i == 1
            % something
        elseif i == Nx
            for i_branch = 1:N_branches % no-flow boundary condition
                u_branch_new(i_branch, Nx_branch) = u_branch_new(i_branch, Nx_branch - 1);
            end
        else
            for ix = 2:(Nx-1) 
                term_excite = magnify * u(ix) * (1 - u(ix)) * (u(ix) - threshold);
                term_2 = 0;
                for i_branch = 1:N_branches
                    term_2 = term_2 + D * (u(ix)-u_branch(i_branch,1)) / Dx^2;
                end
                u_new(ix) = u(ix) + Dt * (term_2 + term_excite);
            end
            for ix = 1:Nx
                v_new(ix) = v(ix) + Dt*(u(ix)-v(ix));
            end
        end
    end
    
    % Enforce Neumann boundary conditions (du/dx=0) on the ends
    % of the system:
    u_new(1) = u_new(2);
    u_new(Nx) = u_new(Nx-1);
    
    % Update for the next timestep:
    u = u_new;
    v = v_new;
    u_branch = u_branch_new;
    v_branch = v_branch_new;
    
    % for plotting
    x_plot = [x, x_branch];
    u_plot = [u, u_branch(i_branch_plot,:)];
    v_plot = [v, v_branch(i_branch_plot,:)];
    
    % Plot every so often:
    if (mod(it,itplot)==0)
        figure(1);
        plot(x,u,x,v);
        axis([x(1),x(Nx),0,1]); % define the plot axes
        title(sprintf('u vs. x at time %f',it*Dt));
        xlabel('x'); ylabel('u');
        plot(x_plot,u_plot,x_plot,v_plot);
        axis([x_plot(1),x_plot(end),0,1]); % define the plot axes
        title(sprintf('u and v vs. x at time %f',it*Dt));
        xlabel('x'); ylabel('u and v');
        drawnow;
    end
    
end