%%% 09/2016 Pierre-Yves Taunay
%%% This script compares the Lax-Friedrich flux in both explicit form with
%%% an Euler time-stepping algorithm, and in semi-explicit form with a TVD
%%% RK3 time-stepping algorithm. 
%%% We are solving the 1D advection equation, with a starting value
%%% featuring smooth functions and discontinuities

clear all

% Define our flux function
% syms q flux_func(q)
% flux_func(q) = q;

xmin = -1;
xmax = 1;
ncell = 200;

tmin = 0;
tmax = 20;

dx = abs(xmax-xmin)/ncell;
xcell = -1+dx/2:dx:1-dx/2;
xcell = xcell';

u0 = uinit(xcell);

CFL = 0.4;
dt = CFL*dx;

up1 = zeros(size(xcell,1),1);

%% Advection with a simple Euler + LF
u = u0;
urk3 = up1;
t = 0;

while t<=tmax
    for i = 1:ncell
        if i == 1
            up1(i,1) = 1/2*(u(ncell,1)+u(i+1,1))-dt/(2*dx)*(u(i+1,1)-u(ncell,1));
        elseif i==ncell
            up1(i,1) = 1/2*(u(i-1,1)+u(1,1))-dt/(2*dx)*(u(1,1)-u(i-1,1));
        else
            up1(i,1) = 1/2*(u(i-1,1)+u(i+1,1))-dt/(2*dx)*(u(i+1,1)-u(i-1,1));
        end
    end
    
    u = up1;
    t = t+dt;
end

plot(xcell,u0,xcell,up1);

% 
% 
%% Advection with RK3 TVD + LF
f1h = zeros(size(xcell,1),2);

u = u0;
t = 0;
while t<=tmax

    u1 = u + dt * lf_lop(f1h,u,ncell,dx);
    u2 = 3/4*u + 1/4*u1 + 1/4*dt*lf_lop(f1h,u1,ncell,dx);
    urk3 = 1/3*u+2/3*u2+2/3*dt*lf_lop(f1h,u2,ncell,dx);
    
    u = urk3;
    t = t+dt;
    
end

hold on
plot(xcell,urk3);

%% EUler + BW
f1h = zeros(size(xcell,1),2);

u = u0;
t = 0;
while t<=tmax

   u1 = u + dt * bw_lop(f1h,u,ncell,dt,dx);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*bw_lop(f1h,u1,ncell,dt,dx);
   urk3 = 1/3*u + 2/3*u2 + 2/3*dt*bw_lop(f1h,u2,ncell,dt,dx);
    up1 = u + dt*bw_lop(f1h,u,ncell,dt,dx);
    u = up1;
    t = t+dt;
    
end

hold on
plot(xcell,u0,xcell,up1);


%% Euler + BW/LF + minmod limiter
f1h = zeros(size(xcell,1),2);

u = u0;
t = 0;
while t<=tmax

    up1 = u + dt * bwlfl_lop(f1h,u,ncell,dt,dx);
    
    u = up1;
    t = t+dt;
    
end

hold on
plot(xcell,up1);

%% WENO 3 FDM
f1h = zeros(size(xcell,1),2);
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno3_lop(u,dx,ncell);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno3_lop(u1,dx,ncell);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno3_lop(u2,dx,ncell);
    
   u = urk3;
        
  
   t = t+dt;
   
    
end

hold on
plot(xcell,u0,'b',xcell,urk3);

%% WENO 5 FDM 
f1h = zeros(size(xcell,1),2);
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_lop(u,dx,ncell);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_lop(u1,dx,ncell);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_lop(u2,dx,ncell);
    
   u = urk3;
        
  
   t = t+dt;
   
    
end


hold on
plot(xcell,u0,'b',xcell,urk3);
% 
%% WENO 3 FVM approach
f1h = zeros(size(xcell,1),2);
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno3_fvm_lop(u,dx,dt,ncell,'FORCE',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno3_fvm_lop(u1,dx,dt,ncell,'FORCE',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno3_fvm_lop(u2,dx,dt,ncell,'FORCE',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

hold on
plot(xcell,u0,xcell,urk3);


%% WENO 5 FVM approach
f1h = zeros(size(xcell,1),2);
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'LW','',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'LW','',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'LW','',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

hold on
plot(xcell,u0,xcell,urk3,'ro');
