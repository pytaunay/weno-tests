function [ un,up ] = weno5Core( u,v )

%% u_{i+1/2}{-}
% Smoothness indicators
vmm = circshift(v,2);
vm  = circshift(v,1);
vp  = circshift(v,-1);
vpp = circshift(v,-2);

% Polynomials
p0 = (2*vmm - 7*vm + 11*v)/6;
p1 = ( -vm  + 5*v  + 2*vp)/6;
p2 = (2*v   + 5*vp - vpp )/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(vmm-2*vm+v  ).^2 + 1/4*(vmm-4*vm+3*v).^2; 
b1 = 13/12*(vm -2*v +vp ).^2 + 1/4*(vm-vp).^2;
b2 = 13/12*(v  -2*vp+vpp).^2 + 1/4*(3*v-4*vp+vpp).^2;

% [b0,b1,b2,p0,p1,p2] = weno5BC(b0,b1,b2,p0,p1,p2,'L',vmm,vm,v,vp,vpp);

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
un = w0.*p0 + w1.*p1 + w2.*p2;

%% % $u_{i-1/2}^{+}$ 
% Choose the negative fluxes, 'u', to compute the left cell boundary flux:
umm = circshift(u,2);
um  = circshift(u,1);
up  = circshift(u,-1);
upp = circshift(u,-2);

% Polynomials
p0 = ( -umm + 5*um + 2*u  )/6;
p1 = ( 2*um + 5*u  - up   )/6;
p2 = (11*u  - 7*up + 2*upp)/6;

% Smooth Indicators (Beta factors)
b0 = 13/12*(umm-2*um+u  ).^2 + 1/4*(umm-4*um+3*u).^2; 
b1 = 13/12*(um -2*u +up ).^2 + 1/4*(um-up).^2;
b2 = 13/12*(u  -2*up+upp).^2 + 1/4*(3*u -4*up+upp).^2;


% [b0,b1,b2,p0,p1,p2] = weno5BC(b0,b1,b2,p0,p1,p2,'R',umm,um,u,up,upp);

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
up = w0.*p0 + w1.*p1 + w2.*p2;
% 
% %% Remove numerical error
% for i = 1:size(up,1)
%     for k = 1:3
%       if(abs(up(i,k)) < 1e-10)
%           up(i,k) = 0;
%       end
%       
%       if(abs(un(i,k)) < 1e-10)
%           un(i,k) = 0;
%       end
% 
% 
% end

