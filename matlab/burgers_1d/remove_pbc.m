function [ unout,upout ] = remove_pbc( un,up,U )

v=U;
u=U;

ncell = size(U,1);
BIG = 1e20;

unout = un;
upout = up;

vmm = circshift(v,2);
vm  = circshift(v,1);
vp  = circshift(v,-1);
vpp = circshift(v,-2);

umm = circshift(u,2);
um  = circshift(u,1);
up  = circshift(u,-1);
upp = circshift(u,-2);


%% Remove pbc on i = 1
% u_{i+1/2}{-}
% Polynomials
p0 = (2*BIG - 7*BIG + 11*v)/6;
p1 = ( -BIG  + 5*v  + 2*vp)/6;
p2 = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(BIG-2*BIG+v  ).^2 + 1/4*(BIG-4*BIG+3*v).^2; 
b1 = 13/12*(BIG -2*v +vp ).^2 + 1/4*(BIG-vp).^2;
b2 = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0 = 1/10; d1 = 6/10; d2 = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
tmpn = w0.*p0 + w1.*p1 + w2.*p2;
unout(1,1) = tmpn(1,1);

% $u_{i-1/2}^{+}$ 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux
% Polynomials
p0 = ( -BIG + 5*BIG + 2*u  )/6;
p1 = ( 2*BIG + 5*u  - up   )/6;
p2 = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(BIG-2*BIG+u  ).^2 + 1/4*(BIG-4*BIG+3*u).^2; 
b1 = 13/12*(BIG -2*u +up ).^2 + 1/4*(BIG-up).^2;
b2 = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0 = 3/10; d1 = 6/10; d2 = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
tmpp = w0.*p0 + w1.*p1 + w2.*p2;

upout(1,1) = tmpp(1,1);

%% Remove PBC at i = 2

% Polynomials
p0 = (2*BIG - 7*vm + 11*v)/6;
p1 = ( -vm  + 5*v  + 2*vp)/6;
p2 = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(BIG-2*vm+v  ).^2 + 1/4*(BIG-4*vm+3*v).^2; 
b1 = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
b2 = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% Constants
d0 = 1/10; d1 = 6/10; d2 = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
tmpn = w0.*p0 + w1.*p1 + w2.*p2;

unout(2,1) = tmpn(2,1);


% $u_{i-1/2}^{+}$ 
% Polynomials
p0 = ( -BIG + 5*um + 2*u  )/6;
p1 = ( 2*um + 5*u  - up   )/6;
p2 = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(BIG-2*um+u  ).^2 + 1/4*(BIG-4*um+3*u).^2; 
b1 = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
b2 = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;

% Constants
d0 = 3/10; d1 = 6/10; d2 = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
tmpp = w0.*p0 + w1.*p1 + w2.*p2;

upout(2,1) = tmpp(2,1);


%% Remove PBC at i = ncell -1

% Polynomials
p0 = (2*vmm - 7*vm + 11*v)/6;
p1 = ( -vm  + 5*v  + 2*vp)/6;
p2 = (2*v   + 5*vp - BIG )/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
b1 = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
b2 = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+BIG).^2;

% Constants
d0 = 1/10; d1 = 6/10; d2 = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
tmpn = w0.*p0 + w1.*p1 + w2.*p2;

unout(ncell-1,1) = tmpn(ncell-1,1);


% $u_{i-1/2}^{+}$ 
% Polynomials
p0 = ( -umm + 5*um + 2*u  )/6;
p1 = ( 2*um + 5*u  - up   )/6;
p2 = (11*u  - 7*up + 2*BIG)/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
b1 = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
b2 = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+BIG).^2;

% Constants
d0 = 3/10; d1 = 6/10; d2 = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
tmpp = w0.*p0 + w1.*p1 + w2.*p2;

upout(ncell-1,1) = tmpp(ncell-1,1);

%% Remove PBC at i = 2

% Polynomials
p0 = (2*vmm - 7*vm + 11*v)/6;
p1 = ( -vm  + 5*v  + 2*BIG)/6;
p2 = (2*v   + 5*BIG - BIG )/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
b1 = 13/12*(vm -2*v +BIG ).^2 + 1/4*(vm-BIG).^2;
b2 = 13/12*(v  -2*BIG+BIG).^2 + 1/4*(3*v-4*BIG+BIG).^2;

% Constants
d0 = 1/10; d1 = 6/10; d2 = 3/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i+1/2}^{-}$;
tmpn = w0.*p0 + w1.*p1 + w2.*p2;

unout(ncell,1) = tmpn(ncell,1);


% $u_{i-1/2}^{+}$ 
% Polynomials
p0 = ( -umm + 5*um + 2*u  )/6;
p1 = ( 2*um + 5*u  - BIG   )/6;
p2 = (11*u  - 7*BIG + 2*BIG)/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
b1 = 13/12*(um -2*u +BIG ).^2 + 1/4*(um-BIG).^2;
b2 = 13/12*(u  -2*BIG+BIG).^2 + 1/4*(3*u -4*BIG+BIG).^2;

% Constants
d0 = 3/10; d1 = 6/10; d2 = 1/10; epsilon = 1e-6;

% Alpha weights 
alpha0 = d0./(epsilon + b0).^2;
alpha1 = d1./(epsilon + b1).^2;
alpha2 = d2./(epsilon + b2).^2;
alphasum = alpha0 + alpha1 + alpha2;

% ENO stencils weigths
w0 = alpha0./alphasum;
w1 = alpha1./alphasum;
w2 = alpha2./alphasum;

% Numerical Flux at cell boundary, $u_{i-1/2}^{+}$;
tmpp = w0.*p0 + w1.*p1 + w2.*p2;

upout(ncell,1) = tmpp(ncell,1);

end

