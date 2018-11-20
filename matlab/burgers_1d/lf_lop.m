function out = lf_lop(f1h,uin,ncell,dx)

    for i=2:ncell-1
        f1h(i,1) = lf_flux(uin(i,1),uin(i+1,1));
        f1h(i,2) = lf_flux(uin(i-1,1),uin(i,1));
    end
    
    i = 1;
    f1h(i,1) = lf_flux(uin(i,1),uin(i+1,1));
    f1h(i,2) = lf_flux(uin(ncell,1),uin(i,1));
    
    i=ncell;
    f1h(i,1) = lf_flux(uin(i,1),uin(1,1));
    f1h(i,2) = lf_flux(uin(i-1,1),uin(i,1));
   
    out = -1/dx*(f1h(:,1)-f1h(:,2));
    
end

