%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Solving 1-D Euler system of equations with 5th order
%          Weighted Essentially Non-Oscilaroty (MOL-WENO5-LF)
%
%        dq_i/dt + df_i/dx = 0, for x \in [a,b] and i =1,. ..,D
%
%           coded by Manuel A. Diaz, manuel.ade'at'gmail.com 
%            Institute of Applied Mechanics, NTU, 2012.08.25
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves the Sod's shock tube problem (IC=1)
%
% t=0                                 t=tEnd
% Density                             Density
%   ****************|                 *********\
%                   |                           \
%                   |                            \
%                   |                             ****|
%                   |                                 |
%                   |                                 ****|
%                   ***************                       ***********
%
% coded by Manuel A. Diaz, 2012.12.27. Last modif: 29.04.2016.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ref: C.-W. Shu, High order weighted essentially non-oscillatory schemes
% for convection dominated problems, SIAM Review, 51:82-126, (2009). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes: 
% 1. A fully conservative finite volume implementation of the method of
% lines (MOL) using WENO5 associated with SSP-RK33 time integration method. 
% 2. Sharpenning of contact discontinuities is NOT implemented here.

%clear; %close all; clc;
global gamma

%% Parameters
CFL     = 0.55;	% CFL number
tFinal	= 0.10;	% Final time
nE      = ncell;  % Number of cells/Elements
gamma   = 1.4;  % Ratio of specific heats for ideal di-atomic gas
IC      = input;	% 10 IC cases are available
plot_fig= 1;

% Discretize spatial domain
a=0; b=1; dx=(b-a)/nE; nx=nE+1; x=linspace(a,b,nx);

% Set IC
[rho0,u0,p0,q0,tFinal,cfl] = Euler_IC1d(x,IC,gamma);
E0 = p0./((gamma-1)*rho0)+0.5*u0.^2;  % Total Energy density
a0 = sqrt(gamma*p0./rho0);            % Speed of sound
q0=[rho0; rho0.*u0; rho0.*E0];        % vec. of conserved properties


% Discretize time domain
lambda0=max(abs(u0)+a0); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

%% Solver Loop
% Load initial condition
q=q0; it=0; dt=dt0; t=0; lambda=lambda0;

%% Solver Loop
while t<tFinal
 
    % RK Initial step
    qo = q;
    
    % 1st stage
    dF=WENO5LF1d(lambda,q,dx);     q = qo-dt*dF; 
    q(:,1)=qo(:,1); q(:,end)=qo(:,end); % Neumann BCs
    
    % 2nd Stage
    dF=WENO5LF1d(lambda,q,dx);     q = 0.75*qo+0.25*(q-dt*dF);
    q(:,1)=qo(:,1); q(:,end)=qo(:,end); % Neumann BCs

    % 3rd stage
    dF=WENO5LF1d(lambda,q,dx);     q = (qo+2*(q-dt*dF))/3;
    q(:,1)=qo(:,1); q(:,end)=qo(:,end); % Neumann BCs
   
    % compute primary properties
    rho=q(1,:); u=q(2,:)./rho; E=q(3,:)./rho; p=(gamma-1)*rho.*(E-0.5*u.^2);
    a=sqrt(gamma*p./rho); if min(p)<0; error('negative pressure found!'); end
    
    % Update dt and time
    lambda=max(abs(u)+a); dt=CFL*dx/lambda; if t+dt>tFinal; dt=tFinal-t; end
    
    % Update time and iteration counter
	t=t+dt; it=it+1;
end

% Calculation of flow parameters
a = sqrt(gamma*p./rho); M = u./a; % Mach number [-]
p_ref = 101325;             % Reference air pressure (N/m^2)
rho_ref= 1.225;             % Reference air density (kg/m^3)
s = 1/(gamma-1)*(log(p/p_ref)+gamma*log(rho_ref./rho)); 
                            % Entropy w.r.t reference condition
ss = log(p./rho.^gamma);    % Dimensionless Entropy
Q = rho.*u;                 % Mass Flow rate per unit area
e = p./((gamma-1)*rho);     % internal Energy
