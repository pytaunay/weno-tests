function q0 = uinit( xcell,problem,GAM )

q0 = zeros(size(xcell,1),3);

if(strcmp(problem,'Sod'))

   
    % rho_L = 1; rho_R = 0.125
    q0(:,1) = 1.*(xcell<=0.5);
    q0(:,1) = q0(:,1) + 0.125*(xcell>0.5);
    
    % Left and right velocities are 0
    
    % Pressure
    P = 1.*(xcell<=0.5); 
    P = P + 0.1*(xcell>0.5);
    
    % And therefore energy
    e0 = 1/(GAM-1)*P./q0(:,1);
    
    q0(:,3) = q0(:,1).*e0;
    
end



end

