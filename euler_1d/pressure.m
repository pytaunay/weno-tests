% Calculates the pressure
function P = pressure(q,GAM,stateOrPhysical)
 
rho = q(:,1);
q2 = q(:,2);
q3 = q(:,3);

if( strcmp(stateOrPhysical,'state') )
    u = q2./rho;
    e0 = q3./rho;
else
    u = q2;
    e0 = q3;
end

P = rho.*(GAM-1).*(e0 - 1/2*u.^2);

end