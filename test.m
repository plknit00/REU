% Funky diffusion code, that solves
% du/dt = D*(d^2u/dx^2) using the forward Euler method.
% MAKING D AN ARRAY NOW RATHER THAN A VARIABLE

% PLOT VALUES
Nx = 300; % Number of grid points in our simulation model
D = ones(1,Nx); % Diffusion coefficient
D_branch = 1.0;
for ix = 50:400  % small region in middle of exit pathway
    if (ix>=123) && (ix <= 129)
        D(ix) = 0.0213;
%     elseif (ix<123) % right exit pathway
%         D(ix) = 0.85;
%     else % left exit pathway
%         D(ix) = 0.801;
    else % rest of exit pathway
         D(ix) = 0.801;
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

N_branches = 9; % number of branches in fan-out used to model atrium
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
v_branch_point_hist2 = nan(1,Nt);
u_branch_point_hist2 = nan(1,Nt);
cell_val = 122; % which cell are we evaluating for history plots
cell_val2 = 123;
v_branch_point_hist3 = nan(1,Nt);
u_branch_point_hist3 = nan(1,Nt);
v_branch_point_hist4 = nan(1,Nt);
u_branch_point_hist4 = nan(1,Nt);
cell_val3 = Nx;
cell_val4 = Nx + 1;
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
        % Should be D(ix - 0.5) but we cant do non-integrer indecies
        left = D(ix) * (u(ix-1) - u(ix)) / Dx^2;
        % Should be D(ix + 0.5) but we cant do non-integrer indecies
        right = D(ix+1) * (u(ix+1) - u(ix)) / Dx^2;
        term_couple = left + right;
        u_new(ix) = u(ix) + Dt*(term_couple + term_excite);
        coupl(ix,it) = term_couple;
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
    coupl(1,it) = term_couple_val;
    
    % only for cell Nx
    u_new(Nx) = u(Nx) + Dt * D(Nx) * (u(Nx-1) - u(Nx))/Dx^2;
    coupl(Nx,it) = D(Nx) * (u(Nx-1) - u(Nx))/Dx^2;
    for i_branch = 1:N_branches
        u_new(Nx) = u_new(Nx) + Dt*D_branch*(u_branch(i_branch,1) - u(Nx))/Dx^2;
        coupl(Nx,it) = coupl(Nx,it) + D_branch*(u_branch(i_branch,1) - u(Nx))/Dx^2;
    end
    threshold2 = (v(Nx) + b(Nx)) / a;
    u_new(Nx) = u_new(Nx) + Dt*magnify*u(Nx)*(1 - u(Nx))*(u(Nx) - threshold2);
    
    % update each branch
    for i_branch = 1:N_branches
        % first cell of every branch
        curr_L = D_branch*(u(Nx)-u_branch(i_branch,1))/Dx^2;
        curr_R = D_branch*(u_branch(i_branch,2)-u_branch(i_branch,1))/Dx^2;
        excite_cell = magnify*u_branch(i_branch,1)*(1-u_branch(i_branch,1))...
            *(u_branch(i_branch,1)-(v_branch(i_branch,1)+b_branch)/a);
        u_branch_new(i_branch,1) = u_branch(i_branch,1)+Dt*(curr_L+curr_R+excite_cell);
        if (i_branch == 1)
            coupl(1 + Nx_branch,it) = curr_L + curr_R;
        end
        % interior cells
        for ix = 2:(Nx_branch - 1)
            term_1 = D_branch*(u_branch(i_branch,ix-1)-2*u_branch(i_branch,ix)+u_branch(i_branch,ix+1))/Dx^2;
            excite = magnify*u_branch(i_branch,ix)*(1-u_branch(i_branch,ix))...
                *(u_branch(i_branch,ix)-(v_branch(i_branch,ix)+b_branch)/a);
            u_branch_new(i_branch,ix)=u_branch(i_branch,ix)+Dt*(term_1+excite);
            if (i_branch == 1)
                coupl(ix + Nx_branch,it) = term_1;
            end
        end
        % last cell of branch
        curr_L = D_branch*(u_branch(i_branch,Nx_branch - 1)-u_branch(i_branch,Nx_branch))/Dx^2;
        excite_cell = magnify*u_branch(i_branch,Nx_branch)*(1-u_branch(i_branch,Nx_branch))...
            *(u_branch(i_branch,Nx_branch)-(v_branch(i_branch,Nx_branch)+b_branch)/a);
        u_branch_new(i_branch,Nx_branch) = u_branch(i_branch,Nx_branch)+Dt*(curr_L+excite_cell);
        for ix = 1:Nx_branch
            v_branch_new(i_branch,ix)=v_branch(i_branch,ix)+Dt*(u_branch(i_branch,ix)-v_branch(i_branch,ix));
        end
    end
    
    % Update for the next timestep:
    u = u_new;
    v = v_new;
    u_branch = u_branch_new;
    v_branch = v_branch_new;
    
    v_branch_point_hist(it) = v(cell_val-3)';
    u_branch_point_hist(it) = u(cell_val-3)';
    v_branch_point_hist2(it) = v(cell_val2+3)';
    u_branch_point_hist2(it) = u(cell_val2+3)';
    cell_val3 = Nx - 10;
    v_branch_point_hist3(it) = v(cell_val3)';
    u_branch_point_hist3(it) = u(cell_val3)';
    cell_val4 = 311;
    v_branch_point_hist4(it) = v_branch(1,10)';
    u_branch_point_hist4(it) = u_branch(1,10)';
    
    v_traces(1:Nx,it) = v;
    v_traces((Nx+1):end,it) = v_branch(1,:)';
    u_traces(1:Nx,it) = u;
    u_traces((Nx+1):end,it) = u_branch(1,:)';
    
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
%         % % Writing each fram to the file
%         % currFrame = getframe(gcf);
%         % writeVideo(vidObj, currFrame);
%     end
end
% close(vidObj);

%% ************ History plots **********************

% figure(2);
% plot((0:(Nt-1))*Dt,u_branch_point_hist,'r','LineWidth',2); hold on;
% plot((0:(Nt-1))*Dt,u_branch_point_hist2,'g','LineWidth',2); hold on;
% %plot((0:(Nt-1))*Dt,v_branch_point_hist,'b','LineWidth',2); hold on;
% %plot((0:(Nt-1))*Dt,v_branch_point_hist2,'k','LineWidth',2);
% %legend('u cell 1', 'v cell 1', 'u cell 2', 'v cell 2');
% legend('u cell 1', 'u cell 2');
% str = sprintf('u & v at cell 1 = %i and cell 2 = %i',cell_val,cell_val2);
% title(str);
% hold off;
% xlabel('Time'); 
% set(gca,'FontSize',16);

% % Branch point
% figure(3);
% plot((0:(Nt-1))*Dt,u_branch_point_hist3,'r','LineWidth',2); hold on;
% plot((0:(Nt-1))*Dt,u_branch_point_hist4,'g','LineWidth',2); hold on;
% %plot((0:(Nt-1))*Dt,v_branch_point_hist3,'b','LineWidth',2); hold on;
% %plot((0:(Nt-1))*Dt,v_branch_point_hist4,'k','LineWidth',2);
% %legend('u cell 1', 'v cell 1', 'u cell 2', 'v cell 2');
% legend('u cell 1', 'u cell 2');
% str = sprintf('u & v at cell 1 = %i and cell 2 = %i',cell_val3,cell_val4);
% title(str);
% hold off;
% xlabel('Time'); 
% set(gca,'FontSize',16);

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
% plot((0:(Nt-1))*Dt,u_branch_point_hist,'b','LineWidth',2);
% hold on;
% plot((0:(Nt-1))*Dt-6.7448,u_branch_point_hist,'g','LineWidth',2);
% xlabel('Time');ylabel('u at branch point');
% str = sprintf('u at cell num = %i, num branches = %i, b = %f',cell_val,N_branches,b(1));
% title(str);
% hold off;

%% ********** Brain Hansen Stuff // Traces **********

figure(5); %Time histories of u(x,t) & v(x,t) vs. t, number of
for ix = 1:8:(Nx+Nx_branch)
%for ix = (Nx-10):(Nx+10)
    plot((0:(Nt-1))*Dt,u_traces(ix,:)-ix*0.05,'b','LineWidth',2); hold on;
    plot((0:(Nt-1))*Dt,v_traces(ix,:)-ix*0.05,'r','LineWidth',2); hold on;
end
hold off;
xlabel('Time','FontSize',20);ylabel('x','FontSize',20);
str = sprintf('u & v for number of branches = %i, b = %f',N_branches,b(1));
title(str,'FontSize',20);

% figure(6); % ****** Time histories of u(x,t) vs. t, number of *****
% for ix = 1:10:(Nx+Nx_branch)
% % for ix = (Nx-20):(Nx+20)
%     plot((0:(Nt-1))*Dt,u_traces(ix,:)-ix*0.05,'b','LineWidth',2); hold on;
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

%% ***************** Velocity graphs *************


% use linear interpolation
time_arr = zeros(1,Nx+Nx_branch-2);
velocity = zeros(1,Nx+Nx_branch-2);
refract = zeros(1,Nx+Nx_branch-2);
time_arr2 = zeros(1,Nx+Nx_branch-2);
velocity2 = zeros(1,Nx+Nx_branch-2);
time_arr3 = zeros(1,Nx+Nx_branch-2);
velocity3 = zeros(1,Nx+Nx_branch-2);
time_arr4 = zeros(1,Nx+Nx_branch-2);
velocity4 = zeros(1,Nx+Nx_branch-2);
% exclude sinus node since it doesn't carry information
for ix = 51:(Nx + Nx_branch - 2)
    lower_bound = 0.5;
    cell1 = ix;
    t_right = 0;
    t_left = 0;
    % plot overlapping stuff
    for it = 1:Nt
        if (u_traces(cell1,it) > lower_bound)
            t_right = it;
            t_left = it - 1;
            break;
        end
        %refractoriness at time of u wave arrival
        % probably in btw 2 times so linearly interpolate
    end
    big_term = (0.5 - u_traces(ix,t_left)) / (u_traces(ix,t_right) - 0.5);
    alpha = big_term / (1 + big_term);
    time_arr(ix) = t_left + alpha;
    % should be velocity (ix + 1/2) but we can only use integer indicies
    % velcoity is change in dist ( = 1) / change in time 
    velocity(ix) = abs(1/(time_arr(ix) - time_arr(ix-1)));
    for it = 27000:Nt
        if (u_traces(cell1,it) > lower_bound)
            t_right = it;
            t_left = it - 1;
            break;
        end
    end
    big_term = (0.5 - u_traces(ix,t_left)) / (u_traces(ix,t_right) - 0.5);
    alpha = big_term / (1 + big_term);
    time_arr2(ix) = t_left + alpha;
    velocity2(ix) = abs(1/(time_arr2(ix) - time_arr2(ix-1)));
    for it = 53000:Nt
        if (u_traces(cell1,it) > lower_bound)
            t_right = it;
            t_left = it - 1;
            break;
        end
    end
    big_term = (0.5 - u_traces(ix,t_left)) / (u_traces(ix,t_right) - 0.5);
    alpha = big_term / (1 + big_term);
    time_arr3(ix) = t_left + alpha;
    velocity3(ix) = abs(1/(time_arr3(ix) - time_arr3(ix-1)));
    for it = 82000:Nt
        if (u_traces(cell1,it) > lower_bound)
            t_right = it;
            t_left = it - 1;
            break;
        end
    end
    big_term = (0.5 - u_traces(ix,t_left)) / (u_traces(ix,t_right) - 0.5);
    alpha = big_term / (1 + big_term);
    time_arr4(ix) = t_left + alpha;
    velocity4(ix) = abs(1/(time_arr4(ix) - time_arr4(ix-1)));
end

figure(8);
plot((0:(Nx+Nx_branch-3))*Dx,velocity(:),'r','LineWidth',2); hold on;
plot((0:(Nx+Nx_branch-3))*Dx,velocity2(:),'b','LineWidth',2); hold on;
plot((0:(Nx+Nx_branch-3))*Dx,velocity3(:),'k','LineWidth',2); hold on;
plot((0:(Nx+Nx_branch-3))*Dx,velocity4(:),'g','LineWidth',2);
str = sprintf('Velocity of Action Potential');
title(str);
hold off;
xlabel('Position'); 
legend('wave 1', 'wave 2', 'wave 3', 'wave 4');
set(gca,'FontSize',16);

% figure(9);
% plot(refract(:),velocity(:),'r*','LineWidth',2);
% str = sprintf('Velocity of Action Potential vs. Refractoriness at the same time');
% title(str);
% hold off;
% xlabel('Refractoriness'); 
% set(gca,'FontSize',16);

