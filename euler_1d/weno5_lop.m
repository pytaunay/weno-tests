function URET = weno5_lop( q,q0,GAM,DX,DT,fluxType,limiterType)

%% Necessary steps before WENO reconstruction
% Roe averages
[qavep,qaven] = roe_average(q,q0,0);

% Local speed of sound at thhe boundaries
ap = speedOfSound(qavep,GAM);
% an = speedOfSound(qavep,GAM);


% Characteristic matrices
[Dp,Rp,Rinvp] = eigMat(qavep,ap,GAM);
% [Dn,Rn,Rinvn] = eigMat(qaven,an,GAM);


% Change of coordinates to characteristics coordinates
rp = changeCoord(Rinvp,q);
% rn = changeCoord(Rinvn,q);

%% WENO reconstruction on characteristics coordinates

% Requires reconstruction of Q at cell boundaries to calculate the flux at
% the boundaries
v = rp; 
u = rp;

[un,up] = weno5Core(u,v);


%% Return to normal coordinates
up = changeCoord(Rp,up);
un = changeCoord(Rp,un);

%% Compute finite volume residual term, df/dx.
% fi+1/2
% Shift for same indexing on both ui+1/2 and ui-1/2
up = circshift(up,-1);

% Shift to obtain the previous value
umn = circshift(un,1);
ump = circshift(up,1);


alphap = eulerMaxEigVal(un,up);
alphan = eulerMaxEigVal(umn,ump);

%URET = numerical_flux( un, up, DT,DX,flux, dflux,fluxType,limiterType ); %fi+1/2
%URET = URET - numerical_flux( umn, ump, DT,DX,flux, dflux,fluxType,limiterType); % fi+1/2 - fi-1/2
    

URET = -1/DX*URET;


end

