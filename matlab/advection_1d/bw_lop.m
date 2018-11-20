function out = bw_lop(f1h,uin,ncell,dt,dx)

    for i=3:ncell
        f1h(i,1) = bw_flux(uin(i-1,1),uin(i,1),dt,dx);
        f1h(i,2) = bw_flux(uin(i-2,1),uin(i-1,1),dt,dx);
    end
    
    i = 1;
    f1h(i,1) = bw_flux(uin(ncell,1),uin(i,1),dt,dx);
    f1h(i,2) = bw_flux(uin(ncell-1,1),uin(ncell,1),dt,dx);
    
    i = 2;
    f1h(i,1) = bw_flux(uin(i-1,1),uin(i,1),dt,dx);
    f1h(i,2) = bw_flux(uin(ncell,1),uin(i-1,1),dt,dx);
       
    out = -1/dx*(f1h(:,1)-f1h(:,2));
    
end

