function [q0,xmin,xmax,xcell,dx,tmax] = uinit( ncell,problem,GAM )

q0 = zeros(ncell,3);

if(strcmp(problem,'Sod'))
    
    xmin = 0;
    xmax = 1;
    tmax = 0.2;
    
    dx = abs(xmax-xmin)/ncell;
    xcell = xmin+dx/2:dx:xmax-dx/2;
    xcell = xcell';

    
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

elseif(strcmp(problem,'Lax'))
    xmin = -5;
    xmax = 5;
    tmax = 0.5;
    
    dx = abs(xmax-xmin)/ncell;
    xcell = xmin+dx/2:dx:xmax-dx/2;
    xcell = xcell';

    
    % rho_L = 0.445, rho_R = 0.5
    q0(:,1) = 0.445.*(xcell<=0);
    q0(:,1) = q0(:,1) + 0.5*(xcell>0);
    
    % u_L = 0.698, u_R = 0
    q0(:,2) = 0.698*(xcell<=0);
    q0(:,2) = q0(:,2) + 0*(xcell>0);
    
    q0(:,2) = q0(:,2).*q0(:,1);
    
    % P_L = 3.528, P_R = 0.571
    P = 3.528*(xcell<=0);
    P = P + 0.571*(xcell>0);
 
    % And therefore energy
    q0(:,3) = 1/(GAM-1)*P + 1/2*q0(:,2).^2./q0(:,1);
      
end



end

