% Calculates the speed of sound
function a = speedOfSound(q,GAM,stateOrPhysical)

rho = q(:,1);

P = pressure(q,GAM,stateOrPhysical);

a = sqrt(GAM.*P./rho);

end