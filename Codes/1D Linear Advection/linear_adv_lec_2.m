% explicit linear advection code 1-D

x = linspace(0,1,101); % create a computational grid between 0 and 1 with 100 grid points
% alternatively you can do this: dx = 0.01; x = 0:dx:100
dx = x(2)-x(1);
Ni = length(x);

u0 = 1; % advection speed;
sigma = 0.05;
Q = exp(-(x-0.5).^2/sigma^2); % initial profile of Q
Q_init = Q; % save the initial profile for the final plot

% impose boundary conditions here - periodic
Q(1) = Q(Ni-1);
Q(Ni) = Q(2);

dt = 0.01; % FIXME: time step should follow wave propagation
Time = 0; % keep track on the time step

figure('position',[442   668   988   280]) % create a blank figure to show the advection results

for n=1:200 % the time stepping loop - let's try advect 200 steps maximum

    Time = Time + dt; % the current time
    
    for i=2:Ni-1 % index stops at Ni-1 because we need i+1
        Q(i) = Q(i) + u0*dt/dx.*(Q(i+1)-Q(i)); % our first forward-difference scheme!
    end
    
    % This is called a periodic boundary condition
    Q(1) = Q(Ni-1);
    Q(Ni) = Q(2);
    
    % plot the Q profile as a function of x
    plot(x,Q);
    xlabel('x'),ylabel('Q'),title(['Simulation Time = ',num2str(Time)]);
    set(gca,'fontsize',14); % make the font better for visualization
    pause(0.1) 
    
    if(Time > 1)
        break; % stop the time-loop if simulation time reaches 1
    end
    
end