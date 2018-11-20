%%% 09/2016 Pierre-Yves Taunay
%%% Function: Jacobian
%%% Input: the vector of physical quantities q and the value of gamma 
%%% for the gas of interest
%%% The function calculates the Jacobian of the flux function based on the
%%% input vector q. 
%%% q fed to the Jacobian function is the result of the Roe averages.
%%% Therefore, it does not requires to modify q to read physical quantities

%%% Jacobian matrix
% System:
% d/dt U + d/dx F(U) = 0 with U = [rho, rho*U, rho*E] = [x1 , x2 , x3]
%           _              _
%          |  rho*u         |
% F(U) =   | rho*u^2 + P    |
%          |rho*u*(E+P/rho) |
%          --             --
% with P = (gamma - 1)*(E - 1/2*rho*u^2)
function J = jacobian( q, GAM )

rho = q(:,1); % rho
u = q(:,2)./rho; % u
e0 = q(:,3)./rho; % e

J = cell(size(q,1),1);

% Populate the array holding all Jacobian matrices
for i = 1:size(q,1)
    Jin = zeros(3,3);
    Jin(1,2) = 1;

    Jin(2,1) = (GAM-3)/2*u(i).^2;
    Jin(2,2) = (3-GAM).*u(i);
    Jin(2,3) = GAM-1;

    Jin(3,1) = (GAM-1)*u(i).^3 - e0(i).*u(i)*GAM;
    Jin(3,2) = e0(i)*GAM + (1-GAM)*3/2*u(i).^2;
    Jin(3,3) = GAM*u(i);
    
    J{i} = Jin;
end
    
end

