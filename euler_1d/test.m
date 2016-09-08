%%% 09/2016
%%% Pierre-Yves Taunay
%%% Test script for Euler equations solved with WENO 5

%% Initial conditions
clear all

problem = 'Sod';
GAM = 1.4;

xmin = 0;
xmax = 1;
ncell = 10;

tmin = 0;
tmax = 0.2;

dx = abs(xmax-xmin)/ncell;
xcell = xmin+dx/2:dx:xmax-dx/2;
xcell = xcell';

q0 = uinit(xcell,problem,GAM);

CFL = 0.4;
dt = CFL*dx;


%% Computations

t = 0;

q = q0;

while ( t < tmax)
    
    % For each RK step
    q1 = q + dt * weno5_lop(q,dx,dt,'LF','');
    q2 = 3/4*q + 1/4*q1 + 1/4*dt*weno5_lop(q1,dx,dt,'LF','');
    qrk3 = 1/3*q + 2/3*q2 + 2/3*dt*weno5_lop(q2,dx,dt,'LF','');
    
    q = qrk3;
    t = t+dt;
    
        % Roe averages

        % Calculate the corresponding matrices

        % Transform to characteristic coordinates

        % WENO on characteristics

        % Transform to local coordinates

        % Apply Riemann solver

end



% %% Lax-Wendroff flux
% u = u0;
% t = 0;
% 
% h = waitbar(0,'Running...');
% while t<=tmax
% 
%    u1 = u + dt * weno5_fvm_lop(u,dx,dt,'LW','',@flux_func,@dflux_func);
%    u2 = 3/4*u + 1/4*u1 + 1/4*dt*weno5_fvm_lop(u1,dx,dt,'LW','',@flux_func,@dflux_func);
%    urk3 = 1/3*u+2/3*u2+2/3*dt*weno5_fvm_lop(u2,dx,dt,'LW','',@flux_func,@dflux_func);
%     
%    u = urk3;
% 
%    t = t+dt;
%  
%        step = ceil(t/dt);
%     waitbar(step/nstep);
%    
% end
% close(h);


ulw = urk3;
