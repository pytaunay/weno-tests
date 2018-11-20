%%% 09/2016
%%% Testing the method order
clear all

xmin = -1;
xmax = 1;

tmin = 0;
tmax = 20;
CFL = 0.4;

nvec = 200:100:1000;
L1 = zeros(size(nvec'));
Linf = zeros(size(L1));

%% Iterations
idx = 1;
for ncell = 200:100:1000
%ncell = 200;

    ncell

    dx = abs(xmax-xmin)/ncell;
    xcell = -1+dx/2:dx:1-dx/2;
    xcell = xcell';

    u0 = uinit(xcell);
    dt = CFL*dx;
    urk3 = zeros(size(xcell,1),1);

    nstep = tmax/dt;
    
    %%% Lax-Friedrich flux
    u = u0;
    t = 0;

    while t<=tmax

       u1 = u + dt * weno5_fvm_lop(u,dx,dt,'LF','',@flux_func,@dflux_func);
       u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'LF','',@flux_func,@dflux_func);
       urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'LF','',@flux_func,@dflux_func);

       u = urk3;

       t = t+dt;

    end

    %%% Calculate the L1 and Linf norms
    L1(idx) = 1/size(urk3,1)*sum( abs(urk3 - u0));
    Linf(idx) = max( abs(urk3 - u0));
    
    idx = idx+1;
end
    