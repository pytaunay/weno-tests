function out = bwlfl_lop( f1h, uin, ncell, dt, dx )

    % Temporary fluxes
    f1hlf = zeros(size(f1h));
    f1hbw = zeros(size(f1h));
    r = zeros(size(f1h));
    
    for i=2:ncell-1
        f1hlf(i,1) = lf_flux(uin(i,1),uin(i+1,1));
        f1hlf(i,2) = lf_flux(uin(i-1,1),uin(i,1));
    end
    
    i = 1;
    f1hlf(i,1) = lf_flux(uin(i,1),uin(i+1,1));
    f1hlf(i,2) = lf_flux(uin(ncell,1),uin(i,1));
    
    i=ncell;
    f1hlf(i,1) = lf_flux(uin(i,1),uin(1,1));
    f1hlf(i,2) = lf_flux(uin(i-1,1),uin(i,1));

    for i=3:ncell
        f1hbw(i,1) = bw_flux(uin(i-1,1),uin(i,1),dt,dx);
        f1hbw(i,2) = bw_flux(uin(i-2,1),uin(i-1,1),dt,dx);
    end
    
    i = 1;
    f1hbw(i,1) = bw_flux(uin(ncell,1),uin(i,1),dt,dx);
    f1hbw(i,2) = bw_flux(uin(ncell-1,1),uin(ncell,1),dt,dx);
    
    i = 2;
    f1hbw(i,1) = bw_flux(uin(i-1,1),uin(i,1),dt,dx);
    f1hbw(i,2) = bw_flux(uin(ncell,1),uin(i-1,1),dt,dx);
    
    % van Leer limiter
    for i=3:ncell-1
        r(i,1) = (uin(i,1) - uin(i-1,1))/(uin(i+1,1)-uin(i,1));
        r(i,2) = (uin(i-1,1) - uin(i-2,1))/(uin(i,1)-uin(i-1,1));
    end
    
    i = 1;
    r(i,1) = (uin(i,1) - uin(ncell,1))/(uin(i+1,1)-uin(i,1));
    r(i,2) = (uin(ncell,1) - uin(ncell-1,1))/(uin(i,1)-uin(ncell,1));
 
    i = 2;
    r(i,1) = (uin(i,1) - uin(i-1,1))/(uin(i+1,1)-uin(i,1));
    r(i,2) = (uin(i-1,1) - uin(ncell,1))/(uin(i,1)-uin(i-1,1));
    
    i = ncell;
    r(i,1) = (uin(i,1) - uin(i-1,1))/(uin(1,1)-uin(i,1));
    r(i,2) = (uin(i-1,1) - uin(i-2,1))/(uin(i,1)-uin(i-1,1));
    
    % minmod
    for k = 1:2
         C = max(0,max(0,min(1,r(:,k)))); % minmod
        % C = max(max(0,min(1,2*r(:,k))),min(r(:,k),2)); superbee
        % C = max(0,min(2*r(:,k), min( (1+r(:,k))/2 , 2))); % MC
        
        f1h(:,k) = f1hlf(:,k) + C.*(f1hbw(:,k)-f1hlf(:,k));
 
  
    end
    
    out = -1/dx*(f1h(:,1)-f1h(:,2));

    
end

