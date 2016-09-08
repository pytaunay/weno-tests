% Calculates the speed of sound
function a = speedOfSound(q,GAM)

rho = q(:,1);

P = pressure(q,GAM);

a = sqrt(GAM.*P./rho);

end