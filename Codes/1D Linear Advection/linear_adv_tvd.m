%%  explicit linear advection code 1-D

x = linspace(0,1,101); % create a computational grid between 0 and 1 with 100 grid points
% alternatively you can do this: dx = 0.01; x = 0:dx:100
dx = x(2)-x(1);
Ni = length(x);

u0 = 1; % advection speed;
sigma = 0.1;
Q = exp(-(x-0.5).^2/sigma^2); % initial profile of Q

Q = x*0;
Q(abs(x-0.5)<0.1)=1;
Q_init = Q; % save the initial profile for the final plot

% impose boundary conditions here - periodic
Q(1) = Q(Ni-1);
Q(Ni) = Q(2);

% dt = 0.001; % FIXME: time step should follow wave propagation
Time = 0; % keep track on the time step

% figure('position',[442   668   988   280]) % create a blank figure to show the advection results
count = 1;
tt = 0;
mm = 0;
Qc = Q;
CFL = 0.8;
dt = abs(dx/u0*CFL);
for n=1:500 % the time stepping loop - let's try advect 200 steps maximum

    % plot the Q profile as a function of x
    if(mod(n,5)==-1) % only plot every 5 steps
        plot(x,Q,'-o');
        hold on
        plot(x,Qc,'-o');
        hold on        
        plot(x,x*0+1,'k--');
        
%         Q_ana = x*0;
%         Q_ana(abs(x+u0*Time-0.5)<=0.11)=1;        
%         plot(x,Q_ana,'r-');
        hold off     
        legend('Backward FD','central-FD') 
        
        xlabel('x'),ylabel('Q'),title(['dt =',num2str(dt),', Simulation Time = ',num2str(Time)]);
        set(gca,'fontsize',14); % make the font better for visualization
        ylim([-0.2 1.5])
        
        
        pause(0.1)
        
%         tt(count) = Time;
%         mm(count) = max(Q);
        % save figure
        str = sprintf('%0*d',4,count);
        image_fn = ['~/Data/solver/png/lec2-',str];
        disp(image_fn);
%         saveas(gcf,image_fn,'png');
        count = count+1;
    end
   
    Time = Time + dt; % the current time
%     
%     for i=2:Ni-1 % index stops at Ni-1 because we need i+1
%         Q(i) = Q(i) + u0*dt/dx.*(Q(i+1)-Q(i)); % our first forward-difference scheme!
%     end
%     
%     for i=2:Ni-1 % index stops at Ni-1 because we need i+1
%         Qc(i) = Qc(i) + u0*dt/(2*dx).*(Qc(i+1)-Qc(i-1)); % central-difference scheme!
%     end

%     for i=2:Ni-1 % index stops at Ni-1 because we need i+1
%         Qc(i) = Qc(i) + u0*dt/(2*dx).*(Qc(i+1)-Qc(i-1)); % lax wendroff scheme!
%     end

% try slope limiting methods
    % first calculate the slope
%     s = (Q(3:end)-Q(1:end-2))/dx/2; % Fromm's
    s = (Q(2:end)-Q(1:end-1))/dx; % beam warming
%     s = (Q(2:end)-Q(1:end-1))/dx; % Lax Wend
%     s(:) = 0;
%     % periodic for slope
%     s(1) = s(Ni-3);
%     s(2) = s(Ni-2);
%     
%     s(Ni-1) = s(3); 
%     s(Ni) = s(4);     
%   

%     for i=2:Ni-1 % index stops at Ni-1 because we need i+1
% %         s(i) = (Q(i+1) - Q(i))/dx;   
% %         Q(i) = Q(i) - u0*dt/dx/2.*(Q(i+1)-Q(i-1)) + 1/2*(u0*dt/dx)^2*(Q(i-1)-2*Q(i)+Q(i+1));
%         
%         
%     end    
    for i=2:Ni-1 % index stops at Ni-1 because we need i+1
%         Q(i) = Q(i) - u0*dt/dx.*(Q(i)-Q(i-1)) - 1/2*u0*dt/dx*(dx -u0*dt)*(s(i) - s(i-1));   
        Q(i) = Q(i) - u0*dt/dx/2.*(Q(i+1)-Q(i-1)) + 1/2*(u0*dt/dx)^2*(Q(i-1)-2*Q(i)+Q(i+1));
        
        
    end    

    % This is called a periodic boundary condition
    Q(1) = Q(Ni-1);
    Q(Ni) = Q(2);    
    
    subplot(1,2,1)
    plot(Q)
    subplot(1,2,2)
    plot(s)
    pause(0.1)
    
%     % diagnostics
%     if(mod(n,5)==0) % only store max Q every 5 steps
%         
%         tt(count) = Time;
%         mm(count) = max(Q);
%         count = count+1;   
%         
%     end
    
    if(Time > 1.)
        break; % stop the time-loop if simulation time reaches 1
    end
    
end
%%
% % plot the initial and final profiles of Q
% figure('position',[442   668   988   280]) % create a blank figure to show the advection results
% plot(x,Q_init,'-ro'); hold on
% plot(x,Q,'-bo')
% Q_ana = x*0;
% Q_ana(abs(x-0.3)<0.1)=1;
% % plot(x,Q_ana,'--b')
% xlabel('x')
% ylabel('Q')
% ylim([-0.1 1.2])
% set(gca,'fontsize',14)
%%
% % plot the initial and final profiles of Q
% figure('position',[442   668   988   280],'color','w') % create a blank figure to show the advection results
% plot(x,Q_init,'-ro'); hold on
% plot(x,Q,'-bo')
% xlabel('x')
% ylabel('Q')
% ylim([-0.1 1.2])
% set(gca,'fontsize',14)
%         legend('Advected Q','Initial Q') 
% 
