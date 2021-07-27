% Test file to try stuff without breaking good code

% PLOT VALUES
Nx = 600; % Number of grid points in our simulation model
D = ones(1,Nx); % Diffusion coefficient
first_cell = 171;
last_cell = 175;
for ix = 50:600  %  exit pathway
    if (ix>=first_cell) && (ix <= last_cell)
        D(ix) = 0.025;
        %D(ix) = 1.0;
    elseif (ix >= 300) % branches
        D(ix) = 9.0;
    else % rest of exit pathway
         D(ix) = 0.806;
    end
end
Dx = 0.1; % Spacing between grid points in our model, larger # makes larger system length
Dt = 0.1*(Dx*Dx)/max(D); % Timestep size (choose to be numerically stable)
Nt = 100000; % Number of timesteps to run
itplot = 100; % Plot every itplot timesteps
x = (0:(Nx-1))*Dx; % Define the coordinates of the gridpoints on the spatial grid
u = 0.01*rand(1,Nx);

% PARAMETERS and ARRAY ALLOCATION
u_new = zeros(1,Nx);
v_new = zeros(1,Nx);
epsilon = 1/50;
magnify = 1/epsilon; % makes firing term large enough
a = 0.8;
b_val = 0.05;
v = zeros(1, Nx); 
b = zeros(1,Nx); % define b for each cell 
for ix = 1:Nx 
    % if we change this value to be smaller, we may not have enough mass to
    % generate a pulse and push it out of sinus node region
    if (ix<50) % if the cell is in the left region, composed of 49 cells
        b(ix) = -0.25; % set the value of b in this cell to avoid neg number
    else % right region
        b(ix) = b_val; % otherwise, set b in this cell to previous value 0.05
    end
end
v_traces = nan(Nx,Nt);
u_traces = nan(Nx,Nt);

% Timestep loop:
for it = 1:Nt
    for ix = 2:(Nx-1) % for all the interior points of the grid...
        % Advance the value of u to the next timestep, u_new, using
        % the forward Euler method in time to represent du/dt,
        % and centered differencing in space, to represent d^2/dx^2:
        threshold = (v(ix) + b(ix)) / a; % cell fires when u > threshold
        term_excite = magnify * u(ix) * (1 - u(ix)) * (u(ix) - threshold);
        % Should be D(ix - 0.5) but we cant do non-integrer indecies
        left = D(ix) * (u(ix-1) - u(ix)) / Dx^2;
        % Should be D(ix + 0.5) but we cant do non-integrer indecies
        right = D(ix+1) * (u(ix+1) - u(ix)) / Dx^2;
        term_couple = left + right;
        u_new(ix) = u(ix) + Dt*(term_couple + term_excite);
    end
    
    for ix = 1:Nx
        v_new(ix) = v(ix) + Dt*(u(ix)-v(ix));
    end
    
    % Enforce Neumann boundary conditions (du/dx=0) on the ends
    % of the system:
    thresh = (v(1) + b(1)) / a; 
    texcite = magnify * u(1) * (1 - u(1)) * (u(1) - thresh);
    right_current = D(2) * (u(2) - u(1)) / Dx^2;
    term_couple_val = right_current;
    u_new(1) = u(1) + Dt*(term_couple_val + texcite);
    
    % only for cell Nx
    u_new(Nx) = u(Nx) + Dt * D(Nx) * (u(Nx-1) - u(Nx))/Dx^2;
    threshold2 = (v(Nx) + b(Nx)) / a;
    u_new(Nx) = u_new(Nx) + Dt*magnify*u(Nx)*(1 - u(Nx))*(u(Nx) - threshold2);
    
    % Update for the next timestep:
    u = u_new;
    v = v_new;
    
    v_traces(1:Nx,it) = v;
    u_traces(1:Nx,it) = u;
    
%     % Plot every so often:
%     if (mod(it,itplot)==0)
%         figure(1);
%         x_plot = [x, x_branch];
%         u_plot = [u, u_branch(i_branch_plot,:)];
%         v_plot = [v, v_branch(i_branch_plot,:)];
%         plot(x_plot,u_plot,'b','LineWidth',2); hold on;
%         plot(x_plot,v_plot,'r','LineWidth',2); hold off;
%         axis([x_plot(1),x_plot(end),0,1]); % define the plot axes
%         title(sprintf('u and v vs. x at time %f',it*Dt));
%         xlabel('x'); ylabel('u and v');
%         drawnow;
%     end
end

%% ******************* Traces *******************

figure(6); %Time histories of u(x,t) & v(x,t) vs. t, number of
for ix = 1:10:Nx
    plot((0:(Nt-1))*Dt,u_traces(ix,:)-ix*0.05,'b','LineWidth',2); hold on;
    plot((0:(Nt-1))*Dt,v_traces(ix,:)-ix*0.05,'r','LineWidth',2); hold on;
end
hold off;
xlabel('Time','FontSize',20);ylabel('x','FontSize',20);
str = sprintf('u & v for b = %f',b(1));
title(str,'FontSize',20);
set(gca,'FontSize',16);

