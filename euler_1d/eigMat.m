function [D,R,Rinv] = eigMat( q, a, GAM )

rho = q(:,1); % rho
u = q(:,2)./rho; % u
e0 = q(:,3)./rho; % e

D = cell(size(q,1),1);
R = cell(size(q,1),1);
Rinv = cell(size(q,1),1);

% Populate the array holding all eigen matrices
for i = 1:size(q,1)
    
    ai = a(i,1);
    ui = u(i,1);
    
    % Eigenvalues
    Di = zeros(3,3);
        
    Di(1,1) = ui;  
    Di(2,2) = ui-ai;
    Di(3,3) = ui+ai;
    
    D{i} = Di;
    
    % Eigenvectors
    Ri = zeros(3,3);
    
    Ri(:,1) = [1,ui,ui^2/2]';
    Ri(:,2) = [1,ui-ai,ui^2/2-ai*ui+1/(GAM-1)*ai^2]';
    Ri(:,3) = [1,ui+ai,ui^2/2+ai*ui+1/(GAM-1)*ai^2]';
    
    R{i} = Ri;
    Rinv{i} = inv(Ri);

    
end



end

