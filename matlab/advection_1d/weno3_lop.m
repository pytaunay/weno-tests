function URET = weno3_lop( U,DX,NCELL )
% Requires reconstruction of Q at cell boundaries to calculate the flux at
% the boundaries
URET = zeros(NCELL,1);
% 
for i = 1:NCELL
        URET(i,1) = WENO3(U,i,NCELL) - WENO3(U,i-1,NCELL);
       URET(i,1) = -URET(i,1)/DX; 
end


%%% Flux splitting
v = 1/2*(U + 1*U); % 1/2 * (flux(U) + df/du*U)
u = 1/2*(U - 1*U);
u = circshift(u,[0,-1]);


%%% u{I+1/2}{R}
% Smoothness indicators
vm  = circshift(v,[0 1]);
vp  = circshift(v,[0 -1]);

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
hn = o0n .* p0n + o1n .* p1n;
 
% %%% u{I+1/2}{L}
% % Choose the negative fluxes, 'u', to compute the left cell boundary flux:
um  = circshift(u,[0 1]);
up  = circshift(u,[0 -1]);

b0p = (u - um).^2;
b1p = (up-u).^2;

p0p = 1/2*vm + 1/2*v;
p1p = 3/2*v - 1/2*vp;

d0 = 2/3; d1 = 1/3;

a0p = d0./(1e-6 + b0p);
a1p = d1./(1e-6 + b1p);

o0p = a0p ./ (a0p + a1p);
o1p = a1p ./ (a0p + a1p);

% Reconstructed flux
hp = o0p .* p0p + o1p .* p1p;


%%% Final flux reconstructed
URET = (hp-circshift(hp,[0,1])+hn-circshift(hn,[0,1]));
URET = -1/DX*URET;

end

