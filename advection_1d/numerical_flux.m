%%% 09/2016 Pierre-Yves Taunay
%%% Function: numerical_flux
%%% Computers the numerical flux of user's choosing
%%% Inputs:
%%% - un, up: reconstructed values at the boundary of interest (un = u_{-}
%%% or u_{L})
%%% - DT, DX: self-explanatory
%%% - fluxFunc, dfluxFunc: a function handle to the analytical flux
%%% function and its derivative, respectively
%%% - fluxType: string to choose a given type of numerical flux
%%% - limiterType: string to choose a given type of limiter for high-order
%%% fluxes

function out = numerical_flux( un, up, DT,DX,fluxFunc, dfluxFunc,fluxType,limiterType )
    
if(strcmp(fluxType,'LF'))
    out = lf_flux(un,up,fluxFunc,dfluxFunc); % fi+1/2
   
elseif( strcmp(fluxType,'FORCE'))
    
    out = force_flux(un,up,DT,DX,fluxFunc);
       
elseif( strcmp(fluxType,'LW'))
    
    out = lw_flux(un, up, fluxFunc, dfluxFunc);
    
elseif( strcmp(fluxType,'FLIC') )
    
    out = flic_flux(un,up,DT,DX,fluxFunc,dfluxFunc,limiterType);
    
end

end

%%% Function: lf_flux
%%% This function calculates the flux across a boundary using the
%%% Lax-Friedrich formulation in a semi-discrete form
function out = lf_flux(phi_L,phi_R,flux,dflux)

    out = 1/2*(flux(phi_L) + flux(phi_R));
    
    alpha = max(dflux(phi_L),dflux(phi_R));
    
    out = out - alpha/2*(phi_R-phi_L);

    
end


%%% Function: force_flux
%%% This function calculates the flux across a boundary usinng the FORCE
%%% formulation in a semi-discrete form
function out = force_flux( un, up,DT,DX,flux )

    uh = 1/2*( up + un - DT/DX*(flux(up) - flux(un)));
    out = 1/4*( flux(up) + 2*flux(uh) + flux(un) - DX/DT*(up-un));
end

%%% Function: lw_flux
function out = lw_flux( un, up, flux, dflux )

    uh = 1/2*(un + up) - 1/2*max(dflux(un),dflux(up))*(flux(up) - flux(un));
    
    out = flux(uh);


end

%%% Function: flic_flux
function out = flic_flux( un,up,DT,DX,fluxFunc, dfluxFunc,limiterType )

    fforce = force_flux(un,up,DT,DX,fluxFunc);
    flw = lw_flux(un,up,fluxFunc,dfluxFunc);
    
   %%% Calculate r
   upm = circshift(up,1);
   upp = circshift(up,-1);
   
   unm = circshift(un,1);
   unp = circshift(un,-1);
   
   rn = (upm-unm)./(up-un);
   rp = (upp-unp)./(up-un);
   
   %%% Limiter
   phi = min(limiter(rn,limiterType),limiter(rp,limiterType));


    %%% Flux
    out = fforce + phi.*(flw - fforce);
   
end




