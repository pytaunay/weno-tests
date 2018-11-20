%%% 09/2016
%%%  Euler flux vector

function Fq = flux( q, GAM )

    q1 = q(:,1);
    q2 = q(:,2);
    q3 = q(:,3);
    
    Fq = zeros(size(q,1),3);
    
    Fq(:,1) = q2;
    Fq(:,2) = q2.^2./(2*q1)*(3-GAM) + (GAM-1)*q3;
    Fq(:,3) = (1-GAM)*q2.^3./(2*q1.^2) + GAM*q3.*q2./q1;  

end

