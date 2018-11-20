function out = bw_flux(um1,u,dt,dx)

    out = u + 1/2*(1-dt/dx)*(u-um1);

end