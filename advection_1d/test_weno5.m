%%% 09/2016
%%% Testing a few different fluxes and limiters with a WENO5 reconstruction
%%% on the 1D advection problem
clear all

xmin = -1;
xmax = 1;
ncell = 200;

tmin = 0;
tmax = 2;

dx = abs(xmax-xmin)/ncell;
xcell = -1+dx/2:dx:1-dx/2;
xcell = xcell';

u0 = uinit(xcell);

CFL = 0.4;
dt = CFL*dx;

up1 = zeros(size(xcell,1),1);

nstep = tmax/dt;

%% Lax-Friedrich flux
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'LF','',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'LF','',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'LF','',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

ulf = urk3;

%% Force flux
u = u0;
t = 0;

h = waitbar(0,'Running...');

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'FORCE','',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'FORCE','',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'FORCE','',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
   
    step = ceil(t/dt);
    waitbar(step/nstep);
       
end
close(h);

uforce = urk3;

%% Lax-Wendroff flux
u = u0;
t = 0;

h = waitbar(0,'Running...');
while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'LW','',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'LW','',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'LW','',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
 
       step = ceil(t/dt);
    waitbar(step/nstep);
   
end
close(h);


ulw = urk3;

%% FLIC flux w/ minmod
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'FLIC','minmod',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'FLIC','minmod',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'FLIC','minmod',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

uflicmm = urk3;

%% FLIC flux w/ superbee
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'FLIC','superbee',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'FLIC','superbee',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'FLIC','superbee',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

uflicsb = urk3;


%% FLIC flux w/ van Leer
u = u0;
t = 0;

while t<=tmax

   u1 = u + dt * weno5_fvm_lop(u,dx,dt,'FLIC','koren',@flux_func,@dflux_func);
   u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'FLIC','koren',@flux_func,@dflux_func);
   urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'FLIC','koren',@flux_func,@dflux_func);
    
   u = urk3;

   t = t+dt;
       
end

uflicvl = urk3;

%% Plot
subplot(2,3,1)
plot(xcell,u0,'b');
hold on
plot(xcell,ulf,'ro');

subplot(2,3,2)
plot(xcell,u0,'b');
hold on
plot(xcell,uforce,'ro');

subplot(2,3,3)
plot(xcell,u0,'b');
hold on
plot(xcell,ulw,'ro');

subplot(2,3,4)
plot(xcell,u0,'b');
hold on
plot(xcell,uflicmm,'ro');

subplot(2,3,5)
plot(xcell,u0,'b');
hold on
plot(xcell,uflicsb,'ro');

subplot(2,3,6)
plot(xcell,u0,'b');
hold on
plot(xcell,uflicvl,'ro');




