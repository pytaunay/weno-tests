%%% 09/2016
%%% Pierre-Yves Taunay
%%% Simple 1D simulation of a plane wave propagating in vacuum

%% Initial conditions
eps0 = 8.85e-12;
mu0 = 4*pi*1e-7; 

c = 3e8;

ncell = 200;

xmin = -10;
xmax = 10;

dx = (xmax-xmin)/ncell;

xcell = xmin + dx/2:dx:xmax-dx/2;

% cfl = 0.5;
dt = 1/2 * dx/c;

nsteps = 100;

q0 = zeros(ncell,2);

sourceLoc = 10;

f0 = 250e6; % 700 MHz

%% Computations
t = 0;
q = q0;

abcHold = zeros(2,2);

for i=1:nsteps
    k = k+1;

    
    q(sourceLoc,1) = q(sourceLoc,1) + sin(2*pi*f0*t);
    
    % For each RK step
    q1 = q + dt * weno5_lop(q,q0,dx,dt,'LF','',wenoType);
    q2 = 3/4*q + 1/4*q1 + 1/4*dt*weno5_lop(q1,q0,GAM,dx,dt,'LF','',wenoType);  
    qrk3 = 1/3*q + 2/3*q2 + 2/3*dt*weno5_lop(q2,q0,GAM,dx,dt,'LF','',wenoType);
    
    % Absorbing BC
    % k = 1
    q(1,1) = abcHold(1,1);
    abcHold(1,1) = abcHold(1,2);
    abcHold(1,2) = q(1,1);
 
    % k = N
    q(end,1) = abcHold(2,1);
    abcHold(2,1) = abcHold(2,2);
    abcHold(2,2) = q(end,1);
    
      
    q = qrk3;    
    t = t+dt;
end