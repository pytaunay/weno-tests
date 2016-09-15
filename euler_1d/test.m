%%% 09/2016
%%% Pierre-Yves Taunay
%%% Test script for Euler equations solved with WENO 5

%% Initial conditions
clear all
clf
close all

GAM = 1.4;

ncell = 500;

xmin = 0;
xmax = 1;

dx = abs(xmax-xmin)/ncell;
xcell = xmin+dx/2:dx:xmax-dx/2;
xcell = xcell';

tmin = 0;

input = 10;

[rho0,u0,p0,q0,tmax,cfl] = Euler_IC1d(xcell,input,GAM);

% Change time step
u = q0(:,2)./q0(:,1);
a = speedOfSound(q0,GAM,'state');
lam = max( abs(u) + a);


dt = cfl*dx/lam;
wenoType = 'component';

%% Storage for animations
QALL = zeros(size(q0,1),size(q0,2),1000);
DTALL = zeros(1000,1);

%% Computations
t = 0;
q = q0;

k = 1;
QALL(:,:,k) = q;
DTALL(1,1) = dt;

while ( t < tmax)

    k = k+1;

    % For each RK step
    q1 = q + dt * weno5_lop(q,q0,GAM,dx,dt,'LF','',wenoType);
    
    q1(1,:) = q0(1,:);
    q1(end,:) = q0(end,:);
    
    
    q2 = 3/4*q + 1/4*q1 + 1/4*dt*weno5_lop(q1,q0,GAM,dx,dt,'LF','',wenoType);
    q2(1,:) = q0(1,:);
    q2(end,:) = q0(end,:);    
    
    
    qrk3 = 1/3*q + 2/3*q2 + 2/3*dt*weno5_lop(q2,q0,GAM,dx,dt,'LF','',wenoType);
    
    qrk3(1,:) = q0(1,:);
    qrk3(end,:) = q0(end,:);
    
      
    q = qrk3;
    
    % Change time step
    u = q(:,2)./q(:,1);
    a = speedOfSound(q,GAM,'state');
    lam = max( abs(u) + a);
    
    dt = cfl*dx/lam;
    
    if( t + dt > tmax)
        dt = tmax-t;
    end
    
    
    t = t+dt;
    
        % Store data
    QALL(:,:,k) = q;
    DTALL(k,1) = dt;

end

qpy = q;
q0py = q0;

%% Second solver
clear q0 q
WENO5LF

%% Exact solutions
[xe,rhoe,ue,pe,ee,t,Me,se] = EulerExact(q0py(1,1),q0py(1,2)/q0py(1,1),pressure(q0py(1,:),GAM,'state'),...
        q0py(end,1),q0py(end,2)/q0py(end,1),pressure(q0py(end,:),GAM,'state'),...
        tmax,GAM,xmin,xmax);
    
    
      
%% Plotting

% s1=subplot(2,3,1); plot(x,rho,'or',xe,rhoe,'k'); xlabel('x(m)'); ylabel('Density (kg/m^3)');
% s2=subplot(2,3,2); plot(x,u,'or',xe,ue,'k'); xlabel('x(m)'); ylabel('Velocity (m/s)');
% s3=subplot(2,3,3); plot(x,p,'or',xe,pe,'k'); xlabel('x(m)'); ylabel('Pressure (Pa)');
% s4=subplot(2,3,4); plot(x,ss,'or',xe,se,'k'); xlabel('x(m)'); ylabel('Entropy/R gas');
% s5=subplot(2,3,5); plot(x,M,'or',xe,Me,'k'); xlabel('x(m)'); ylabel('Mach number');
% s6=subplot(2,3,6); plot(x,e,'or',xe,ee,'k'); xlabel('x(m)'); ylabel('Internal Energy (kg/m^2s)');

figure
plot(xe,rhoe,xcell,qpy(:,1),x,rho);    
title('Density')

figure
plot(xe,ue,xcell,qpy(:,2)./qpy(:,1),x,u);
title('Velocity')

figure
plot(xe,pe,xcell,pressure(qpy,GAM,'state'),x,p);
title('Pressure')

figure
plot(xe,ee,xcell,qpy(:,3)./qpy(:,1)-1/2*(qpy(:,2)./qpy(:,1)).^2,x,e);
title('Energy')

% figure
% plot(x,p);
% hold on
% plot(xcell,pressure(q,GAM,'state'),'ro');
% 
% 
% figure
% plot(x,ux)
% hold on
% plot(xcell,q(:,2)./q(:,1),'ro');



