% Calculates the pressure
function P = pressure(q,GAM)
    
rho = q(1);
q2 = q(2);
q3 = q(3);

u = q2./rho;
e0 = q3./rho;

P = rho.*(GAM-1).*(e0 - 1/2*u.^2);

end