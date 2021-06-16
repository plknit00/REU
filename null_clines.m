% Phase-plane trajectories of Barkley model

% Parameters
epsilon = 0.02; % st 1/epsilon = 50
a = 0.8;
b = 0.05;

Dt = 1.e-4; % timestep size
Nt = 100000; % number of timesteps
it_mark = 300; % mark trajectories for every it_mark timesteps

% Plotting u-nullcline
figure(3); clf;
subplot(2,1,1);

v_plot = [-0.2,1.2];
plot([0,0],v_plot,'b'); hold on;
grid;
plot([1,1],v_plot,'b');
u_plot = (v_plot+b) / a;
plot(u_plot,v_plot,'b');

% Plotting v-nullcline
u_plot = v_plot;
plot(u_plot,v_plot,'r');

xlabel('u');
ylabel('v');
axis([-0.2,1.2,-0.2,1.2]);
disp('Click on  a point in the graph to set initial conditions');
disp('Type Control-C to exit');

while (1)
    % read mouse location
    [u0,v0] = ginput(1);
    
    % Solve the Barkley equations (w/o diffusion term)
    % (u0,v0) = initial conditions
    u = nan(1,Nt+1);
    v = nan(1,Nt+1);
    u(1) = u0;
    v(1) = v0;
    for it = 1:Nt
        u_new = u(it) + Dt/epsilon*u(it)*(1-u(it))*(u(it)-(v(it)+b)/a);
        v_new = v(it) + Dt*(u(it)-v(it));
        u(it+1) = u_new;
        v(it+1) = v_new;
    end
    % mark initial conditions w/ * and plot trajectory
    subplot(2,1,1)
    plot(u0,v0,'k*');
    plot(u,v,'k');
    % mark trajectory at points that are equally spaced in time
    plot(u(1:it_mark:end),v(1:it_mark:end),'ko');
    % plot trajectors vs. time
    subplot(2,1,2);
    plot((0:Nt)*Dt,u,'b'); hold on;
    plot((0:Nt)*Dt,v,'r');
    axis([0,Nt*Dt,-0.2,1.2]);
    legend('u(t)','v(t)');
    xlabel('time');
    hold off;
end

