function URET = weno3_fvm_lop( U,DX,DT,NCELL,numFlux,flux,dFlux )
% Requires reconstruction of Q at cell boundaries to calculate the flux at
% the boundaries
URET = zeros(NCELL,1);
% % 
% for i = 1:NCELL
%         URET(i,1) = WENO3(U,i,NCELL) - WENO3(U,i-1,NCELL);
%        URET(i,1) = -URET(i,1)/DX; 
% end
v = U;
u = U;

%%% u{I+1/2}{R}
% Smoothness indicators
vm  = circshift(v,1);
vp  = circshift(v,-1);

b0n = (v - vm).^2;
b1n = (vp-v).^2;

p0n = -1/2*vm + 3/2*v;
p1n = 1/2*v + 1/2*vp;

d0 = 1/3; d1 = 2/3;

a0n = d0./(1e-6 + b0n);
a1n = d1./(1e-6 + b1n);

o0n = a0n ./ (a0n + a1n);
o1n = a1n ./ (a0n + a1n);

% Reconstructed flux
un = o0n .* p0n + o1n .* p1n;
 
% %%% u{I-1/2}{L}
% % Choose the negative fluxes, 'u', to compute the left cell boundary flux:
um  = circshift(u,1);
up  = circshift(u,-1);

b0p = (u - um).^2;
b1p = (up-u).^2;

p0p = 1/2*um + 1/2*u;
p1p = 3/2*u - 1/2*up;

d0 = 2/3; d1 = 1/3;

a0p = d0./(1e-6 + b0p);
a1p = d1./(1e-6 + b1p);

o0p = a0p ./ (a0p + a1p);
o1p = a1p ./ (a0p + a1p);

% Reconstructed flux
up = o0p .* p0p + o1p .* p1p;

%%% Final flux reconstructed
% Shift to have same boundaries on both un and up
up = circshift(up,-1);

% Shift to obtain the previous value
umn = circshift(un,1);
ump = circshift(up,1);


if(strcmp(numFlux,'LF'))
    URET = 1/2*( flux(un)+flux(up) - 1*(up-un)); % fi+1/2
    URET = URET - 1/2*( flux(umn)+flux(ump) - 1*(ump-umn)); % fi+1/2 - fi-1/2
    
elseif( strcmp(numFlux,'FORCE'))
    
    uh = 1/2*( up + un - DT/DX*(flux(up) - flux(un)));
    URET = 1/4*( flux(up) + 2*flux(uh) + flux(un) - DX/DT*(up-un));
    
    uh = 1/2*( ump + umn - DT/DX*(flux(ump) - flux(umn)));
    URET = URET - 1/4*( flux(ump) + 2*flux(uh) + flux(umn) - DX/DT*(ump-umn));
end
    
URET = -1/DX*URET;



end
