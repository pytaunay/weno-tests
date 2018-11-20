function URET = weno5_lop( q,q0,GAM,DX,DT,fluxType,limiterType,compOrCharac)

URET = zeros(size(q,1),3);

if(strcmp(compOrCharac,'characteristic'))

    %%% Necessary steps before WENO reconstruction
    % Roe averages w/o PBC
    % Note that the output is in terms of physical quantities, NOT in terms of
    % state vector quantities
     [qavep,qaven] = roe_average(q,q0,0);

    % Local speed of sound at the boundaries
    ap = speedOfSound(qavep,GAM,'phys');
    an = speedOfSound(qaven,GAM,'phys');

    % Characteristic matrices
    [Dp,Rp,Rinvp] = eigMat(qavep,ap,GAM);
    [Dn,Rn,Rinvn] = eigMat(qaven,an,GAM);

    % Change of coordinates to characteristics coordinates
    rp = changeCoord(Rinvp,q);
    rn = changeCoord(Rinvn,q);

    %%% WENO reconstruction on characteristics coordinates
    % Requires reconstruction of Q at cell boundaries to calculate the flux at
    % the boundaries
    % Reconstruct v_{i+1/2}{+/-}
    v = rp; 
    u = rp;

    [vn,vp] = weno5Core(u,v);

    %%% Return to normal coordinates
    up = changeCoord(Rp,vp);
    un = changeCoord(Rp,vn);

    %%% Compute finite volume residual term, df/dx.
    % fi+1/2
    % Shift for same indexing on both ui+1/2 and ui-1/2
     up=circshift(up,-1);
     umn = circshift(un,1);
     ump = circshift(up,1);
    
elseif(strcmp(compOrCharac,'component'))
    %%% Component-wise reconstruction
    v = q;
    u = q;
    [un,up] = weno5Core(u,v);
    
    up = circshift(up,-1);

    % Shift to obtain the previous value
    umn = circshift(un,1);
    ump = circshift(up,1);   

else
    up = circshift(q,-1);
    un = q;
    
    umn = circshift(q,1);
    
    ump = q; 
   
end
    
%%% Enforce boundary conditions
%  up(end,:) = q0(end,:);
%  un(end,:) = q0(end,:);
%  umn(1,:) = q0(1,:);
%  ump(1,:) = q0(1,:);

%% Perform flux calculations
alphap = eulerMaxEigVal(un,up,GAM);
alphan = eulerMaxEigVal(umn,ump,GAM);

Fqp = flux(up,GAM);
Fqn = flux(un,GAM);

Fqmp = flux(ump,GAM);
Fqmn = flux(umn,GAM);

for k = 1:3
    % Rusanov flux
    URET(:,k) = 1/2*( Fqp(:,k) + Fqn(:,k) - alphap.*(up(:,k)-un(:,k)));
    URET(:,k) = URET(:,k) - 1/2*( Fqmp(:,k) + Fqmn(:,k) - alphan.*(ump(:,k)-umn(:,k)));
end


URET = -1/DX*URET;


end

