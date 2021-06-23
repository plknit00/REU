% Standard diffusion code, that solves
% du/dt = D*(d^2u/dx^2) using the forward Euler method.

% PLOT VALUES
D = 1.0; % Diffusion coefficient
Nx = 300; % Number of grid points in our simulation model
Dx = 0.1; % Spacing between grid points in our model, larger # makes larger system length
Dt = 0.1*(Dx*Dx)/D; % Timestep size (choose to be numerically stable)
Nt = 60000; % Number of timesteps to run
itplot = 100; % Plot every itplot timesteps
x = (0:(Nx-1))*Dx; % Define the coordinates of the gridpoints on the spatial grid
u = 0.01*rand(1,Nx);
cell_val = 0; % which cell are we evaluating for history plots

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
        b(ix) = -0.125; % set the value of b in this cell to avoid neg number
    else % right region
        b(ix) = b_val; % otherwise, set b in this cell to previous value 0.05
    end
end

N_branches = 10; % number of branches in fan-out used to model atrium
Nx_branch = 300; % number of cells in each branch
i_branch_plot = 1; % which branch to plot
b_branch = 0.05; % b_value for branches
x_branch = (Nx:(Nx + Nx_branch - 1)) * Dx; % extend coordinate system
% define arrays to keep track of u and v in each cell of each branch
u_branch = 0.01 * rand(N_branches, Nx_branch);
v_branch = zeros(N_branches, Nx_branch);
u_branch_new = zeros(N_branches, Nx_branch);
v_branch_new = zeros(N_branches, Nx_branch);
% history plots 
v_branch_point_hist = nan(1,Nt);
u_branch_point_hist = nan(1,Nt);
% for Brian hansen traces
v_traces = nan(Nx_branch+Nx,Nt);
u_traces = nan(Nx_branch+Nx,Nt);
coupl = zeros(Nx_branch+Nx,Nt); % make a plot

% % To make video
% vidObj = VideoWriter('traces.avi');
% open(vidObj);

% Timestep loop:
for it = 1:Nt
    for ix = 2:(Nx-1) % for all the interior points of the grid...
        % Advance the value of u to the next timestep, u_new, using
        % the forward Euler method in time to represent du/dt,
        % and centered differencing in space, to represent d^2/dx^2:
        threshold = (v(ix) + b(ix)) / a; % cell fires when u > threshold
        term_excite = magnify * u(ix) * (1 - u(ix)) * (u(ix) - threshold);
        term_couple = D * (u(ix-1) - 2*u(ix) + u(ix+1)) / Dx^2;
        u_new(ix) = u(ix) + Dt*(term_couple + term_excite);
        coupl(ix,it) = term_couple;
    end
    for ix = 1:Nx
        v_new(ix) = v(ix) + Dt*(u(ix)-v(ix));
    end
    
    % Enforce Neumann boundary conditions (du/dx=0) on the ends
    % of the system:
    u_new(1) = u_new(2);
    coupl(1,it) = coupl(2,it);
    
    % only for cell Nx
    u_new(Nx) = u(Nx) + Dt * D * (u(Nx-1) - u(Nx))/Dx^2;
    coupl(Nx,it) = D * (u(Nx-1) - u(Nx))/Dx^2;
    for i_branch = 1:N_branches
        u_new(Nx) = u_new(Nx) + Dt*D*(u_branch(i_branch,1) - u(Nx))/Dx^2;
        coupl(Nx,it) = coupl(Nx,it) + D*(u_branch(i_branch,1) - u(Nx))/Dx^2;
    end
    threshold2 = (v(Nx) + b(Nx)) / a;
    u_new(Nx) = u_new(Nx) + Dt*magnify*u(Nx)*(1 - u(Nx))*(u(Nx) - threshold2);
    
    % update each branch
    for i_branch = 1:N_branches
        % first cell of every branch
        curr_L = D*(u(Nx)-u_branch(i_branch,1))/Dx^2;
        curr_R = D*(u_branch(i_branch,2)-u_branch(i_branch,1))/Dx^2;
        excite_cell = magnify*u_branch(i_branch,1)*(1-u_branch(i_branch,1))...
            *(u_branch(i_branch,1)-(v_branch(i_branch,1)+b_branch)/a);
        u_branch_new(i_branch,1) = u_branch(i_branch,1)+Dt*(curr_L+curr_R+excite_cell);
        if (i_branch == 1)
            coupl(1 + Nx_branch,it) = curr_L + curr_R;
        end
        % interior cells
        for ix = 2:(Nx_branch - 1)
            term_1 = D*(u_branch(i_branch,ix-1)-2*u_branch(i_branch,ix)+u_branch(i_branch,ix+1))/Dx^2;
            excite = magnify*u_branch(i_branch,ix)*(1-u_branch(i_branch,ix))...
                *(u_branch(i_branch,ix)-(v_branch(i_branch,ix)+b_branch)/a);
            u_branch_new(i_branch,ix)=u_branch(i_branch,ix)+Dt*(term_1+excite);
            if (i_branch == 1)
                coupl(ix + Nx_branch,it) = term_1;
            end
        end
        % boundary condition
        u_branch_new(i_branch,Nx_branch)=u_branch_new(i_branch,Nx_branch-1);
        for ix = 1:Nx_branch
            v_branch_new(i_branch,ix)=v_branch(i_branch,ix)+Dt*(u_branch(i_branch,ix)-v_branch(i_branch,ix));
        end
    end
    
    % Update for the next timestep:
    u = u_new;
    v = v_new;
    u_branch = u_branch_new;
    v_branch = v_branch_new;
    
    cell_val = Nx - 1;
    v_branch_point_hist(it) = v(cell_val)';
    u_branch_point_hist(it) = u(cell_val)';
    
    v_traces(1:Nx,it) = v;
    v_traces((Nx+1):end,it) = v_branch(1,:)';
    u_traces(1:Nx,it) = u;
    u_traces((Nx+1):end,it) = u_branch(1,:)';
    
    % Plot every so often:
    if (mod(it,itplot)==0)
        figure(1);
        x_plot = [x, x_branch];
        u_plot = [u, u_branch(i_branch_plot,:)];
        v_plot = [v, v_branch(i_branch_plot,:)];
        plot(x_plot,u_plot,'LineWidth',2); hold on;
        plot(x_plot,v_plot,'LineWidth',2); hold off;
        axis([x_plot(1),x_plot(end),0,1]); % define the plot axes
        title(sprintf('u and v vs. x at time %f',it*Dt));
        xlabel('x'); ylabel('u and v');
        drawnow;
        % % Writing each fram to the file
        % currFrame = getframe(gcf);
        % writeVideo(vidObj, currFrame);
    end
end
% close(vidObj);

%% ************ History plots **********************

% figure(2);
% plot((0:(Nt-1))*Dt,v_branch_point_hist,'r','LineWidth',2);
% hold on;
% plot((0:(Nt-1))*Dt,u_branch_point_hist,'b','LineWidth',2);
% str = sprintf('u & v at cell num = %i, num branches = %i, b = %f',cell_val,N_branches,b(1));
% title(str);
% hold off;
% xlabel('Time');
% legend('v', 'u');
% 
% figure(3);
% plot((0:(Nt-1))*Dt,v_branch_point_hist,'r','LineWidth',2);
% hold on;
% plot((0:(Nt-1))*Dt-6.7448,v_branch_point_hist,'g','LineWidth',2);
% xlabel('Time');ylabel('v');
% str = sprintf('v at cell num = %i, num branches = %i, b = %f',cell_val,N_branches,b(1));
% title(str);
% hold off;
% 
% figure(4);
% plot((0:(Nt-1))*Dt,u_branch_point_hist,'r','LineWidth',2);
% hold on;
% plot((0:(Nt-1))*Dt-6.7448,u_branch_point_hist,'g','LineWidth',2);
% xlabel('Time');ylabel('u at branch point');
% str = sprintf('u at cell num = %i, num branches = %i, b = %f',cell_val,N_branches,b(1));
% title(str);
% hold off;

%% ********** Brain Hansen Stuff // Traces **********

% figure(5); %Time histories of u(x,t) & v(x,t) vs. t, number of
% for ix = 1:10:(Nx+Nx_branch)
% %for ix = (Nx-10):(Nx+10)
%     plot((0:(Nt-1))*Dt,u_traces(ix,:)-ix*0.05,'r','LineWidth',2); hold on;
%     plot((0:(Nt-1))*Dt,v_traces(ix,:)-ix*0.05,'b','LineWidth',2); hold on;
% end
% hold off;
% xlabel('Time');ylabel('x');
% str = sprintf('u & v for number of branches = %i, b = %f',N_branches,b(1));
% title(str);
% 
% figure(6); % ****** Time histories of u(x,t) vs. t, number of *****
% for ix = 1:10:(Nx+Nx_branch)
% % for ix = (Nx-20):(Nx+20)
%     plot((0:(Nt-1))*Dt,u_traces(ix,:)-ix*0.05,'r','LineWidth',2); hold on;
% end
% hold off;
% xlabel('Time');ylabel('x');
% str = sprintf('u for number of branches = %i, b = %f',N_branches,b(1));
% title(str);
% 
% figure(7); % ****** Time histories of v(x,t) vs. t, number of *****
% for ix = 1:10:(Nx+Nx_branch)
% % for ix = (Nx-20):(Nx+20)
%     plot((0:(Nt-1))*Dt,v_traces(ix,:)-ix*0.05,'r','LineWidth',2); hold on;
% end
% hold off;
% xlabel('Time');ylabel('x');
% str = sprintf('v for number of branches = %i, b = %f',N_branches,b(1));
% title(str);

%% *********** Overlapping traces plots ************

% figure(8); % ****** Time histories of u(x,t) vs. t, number of *****
% %for ix = 1:10:(Nx+Nx_branch)
% for ix = (Nx-20):10:(Nx+20)
%     plot((0:(Nt-1))*Dt,u_traces(ix,:),'LineWidth',2); hold on;
% end
% hold off;
% xlabel('Time');ylabel('x');
% %set(gca,'FontSize',16);
% str = sprintf('u for number of branches = %i, b = %f',N_branches,b(1));
% title(str);
%  
% figure(9); % ****** Time histories of v(x,t) vs. t, number of *****
% for ix = 1:100:(Nx+Nx_branch)
% %for ix = (Nx-20):(Nx+20)
%     plot((0:(Nt-1))*Dt,v_traces(ix,:),'LineWidth',2); hold on;
% end
% hold off;
% xlabel('Time');ylabel('x');
% str = sprintf('v for number of branches = %i, b = %f',N_branches,b(1));
% title(str);


%% *********** Phase plane plot of specific cell **************
% 
% figure(10);
% % plot nullclines
% cell_num = 100;
% v_phase = [0,1];
% plot([0,0], v_phase, 'b','LineWidth',2); hold on;
% plot([1,1], v_phase, 'b','LineWidth',2);
% if (cell_num <= Nx)
%     b_phase = b(cell_num);
% else 
%     b_phase = b_branch;
% end
% plot((v_phase + b_phase)/a, v_phase, 'b','LineWidth',2);
% plot(v_phase, v_phase, 'r','LineWidth',2);
% xlabel('u'); ylabel('v');
% title(sprintf('Phase plane plot for cell num = %i, num branches  = %i, b = %f', cell_num, N_branches, b(1)));
% comet(u_traces(cell_num,:), v_traces(cell_num,:));
% hold off;

%% ************ Contour Plots // 3D movies **************

% figure(11);
% title(sprintf('u for num branches = %i, b = %f', N_branches, b(1)));
% % vidObj = VideoWriter('traces.avi');
% % open(vidObj);
% for theta = 0:359
%     surf(u_traces(1:10:end,1:1000:end));
%     shading interp;
%     axis vis3d;
%     view(theta, 45);
%     drawnow;
%     % % Writing each fram to the file
%     % currFrame = getframe(gcf);
%     % writeVideo(vidObj, currFrame);
% end
% % close(vidObj);
% 
% figure(12);
% % vidObj = VideoWriter('traces.avi');
% % open(vidObj);
% for theta = 0:359
%     surf(v_traces(1:10:end,1:100:end));
%     title(sprintf('v for num branches = %i, b = %f', N_branches, b(1)));
%     shading interp;
%     axis vis3d;
%     view(theta, 45);
%     drawnow;
%     % % Writing each fram to the file
%     % currFrame = getframe(gcf);
%     % writeVideo(vidObj, currFrame);
% end
% close(vidObj);

%% ************ Stationary Contour Plots ****************

% figure(13);
% contour(u_traces(1:10:end, 1:1000:end))
% colorbar
% title(sprintf('u for num branches = %i, b = %f', N_branches, b(1)));
% 
% figure(14);
% contour(v_traces(1:10:end, 1:1000:end))
% colorbar
% title(sprintf('v for num branches = %i, b = %f', N_branches, b(1)));
% 
% figure(15);
% imagesc(u_traces(1:10:end, 1:1000:end))
% colorbar
% title(sprintf('u for num branches = %i, b = %f', N_branches, b(1)));
% 
% figure(16);
% imagesc(v_traces(1:10:end, 1:1000:end))
% colorbar
% title(sprintf('v for num branches = %i, b = %f', N_branches, b(1)));
% 
% figure(17);
% surf(u_traces(:,1:1000:end))
% shading interp;
% title(sprintf('u for num branches = %i, b = %f', N_branches, b(1)));
% 
% figure(18);
% surf(v_traces(:,1:1000:end))
% shading interp;
% title(sprintf('v for num branches = %i, b = %f', N_branches, b(1)));

%% ************* Time loop of overlapping waves ****************

% figure(19);
% for it = 10301:100:Nt
%     % normal pulses
%     plot((0:((Nx_branch + Nx)-1))*Dx,u_traces(:,it),'r','LineWidth',2); hold on;
%     plot((0:((Nx_branch + Nx)-1))*Dx,v_traces(:,it),'b','LineWidth',2); hold on;
%     % lagging pulses
%     plot((0:((Nx_branch + Nx)-1))*Dx,u_traces(:,it - 10300),'m','LineWidth',2); hold on;
%     plot((0:((Nx_branch + Nx)-1))*Dx,v_traces(:,it - 10300),'c','LineWidth',2); hold off;
%     title(sprintf('u and v vs. x at normal time = %f and lag time = %f',it*Dt,it*Dt-10.30));
%     xlabel('x'); ylabel('u and v');
%     legend('u normal', 'v normal', 'u lag', 'v lag');
%     drawnow;
% end

%% **************** Plot of coupling term ******************

% figure(20);
% % plot of single cell's coupling term in time
% cell_num_val = 250;
% plot((0:(Nt-1))*Dt,coupl(cell_num_val,:),'LineWidth',2);
% xlabel('time'); ylabel('Coupling term value');
% title(sprintf('Coupling term value of cell = %f vs. time',cell_num_val));

figure(21); 
% traces of all cell's coupling terms
for ix = 1:10:(Nx+Nx_branch)
%for ix = (Nx-10):(Nx+10)
    plot((0:(Nt-1))*Dt,coupl(ix,:)-ix*0.25,'r'); hold on;
end
hold off;
xlabel('Time');ylabel('x and coupling term');
title(sprintf('Coupling term value vs. time'));





