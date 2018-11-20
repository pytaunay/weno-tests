function URET = weno5_lop( U,DX,NCELL )
% Requires reconstruction of Q at cell boundaries to calculate the flux at
% the boundaries
URET = zeros(NCELL,1);

%%% Flux splitting
v = 1/2*(U + 1*U); % 1/2 * (flux(U) + df/du*U)
u = 1/2*(U - 1*U);
u = circshift(u,[0,-1]);


%%% Smoothness indicators
vmm = circshift(v,[0 2]);
vm  = circshift(v,[0 1]);
vp  = circshift(v,[0 -1]);
vpp = circshift(v,[0 -2]);

% Polynomials
p0n = (2*vmm - 7*vm + 11*v)/6;
p1n = ( -vm  + 5*v  + 2*vp)/6;
p2n = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
B0n = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
B1n = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
B2n = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0n = 1/10; d1n = 6/10; d2n = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0n = d0n./(epsilon + B0n).^2;
alpha1n = d1n./(epsilon + B1n).^2;
alpha2n = d2n./(epsilon + B2n).^2;
alphasumn = alpha0n + alpha1n + alpha2n;

% ENO stencils weigths
w0n = alpha0n./alphasumn;
w1n = alpha1n./alphasumn;
w2n = alpha2n./alphasumn;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
hn = w0n.*p0n + w1n.*p1n + w2n.*p2n;

%% Left Flux 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
% $u_{i-1/2}^{+}$ 
umm = circshift(u,[0 2]);
um  = circshift(u,[0 1]);
up  = circshift(u,[0 -1]);
upp = circshift(u,[0 -2]);

% Polynomials
p0p = ( -umm + 5*um + 2*u  )/6;
p1p = ( 2*um + 5*u  - up   )/6;
p2p = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
B0p = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
B1p = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
B2p = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0p = 3/10; d1p = 6/10; d2p = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0p = d0p./(epsilon + B0p).^2;
alpha1p = d1p./(epsilon + B1p).^2;
alpha2p = d2p./(epsilon + B2p).^2;
alphasump = alpha0p + alpha1p + alpha2p;

% ENO stencils weigths
w0p = alpha0p./alphasump;
w1p = alpha1p./alphasump;
w2p = alpha2p./alphasump;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
hp = w0p.*p0p + w1p.*p1p + w2p.*p2p;

%% Compute finite volume residual term, df/dx.
URET = (hp-circshift(hp,[0,1])+hn-circshift(hn,[0,1]));
URET = -1/DX*URET;


% % WENO 5th order requires access to i-2, i-1, i, i+1, and i+2
% % Have to be careful with boundary conditions
% fU1 = 0;
% fU2 = 0;
% for i = 1:NCELL
% 
%     fpp = WENO5(U,i,'L',NCELL);
%     fpm = WENO5(U,i+1,'R',NCELL);
%     
%     fmp = WENO5(U,i-1,'L',NCELL);
%     fmm = WENO5(U,i,'R',NCELL);
%     
%    
% %   URET(i,1) = lf_flux(U1,U2);
% %   URET(i,1) = -URET(i,1)/DX;
%  
%     URET(i,1) = fpp + fpm - (fmp+fmm);
%        
% end
% 
% URET = -URET/DX;

end

