%%% 09/2016 Pierre-Yves Taunay
%%% This script compares the Lax-Friedrich flux in both explicit form with
%%% an Euler time-stepping algorithm, and in semi-explicit form with a TVD
%%% RK3 time-stepping algorithm. 
%%% We are solving the 1D advection equation, with a starting value
%%% featuring smooth functions and discontinuities

clear all

xmin = -1;
xmax = 1;
ncell = 80;

tmin = 0;
tmax = 1.1;

dx = abs(xmax-xmin)/ncell;
xcell = -1+dx/2:dx:1-dx/2;
xcell = xcell';

u0 = uinit(xcell);

CFL = 0.8;
dt = CFL*dx;

up1 = zeros(size(xcell,1),1);


%% WENO 5 FVM approach
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'LF','',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'LF','',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'LF','',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

hold on
plot(xcell,u0,xcell,urk3,'ro');
