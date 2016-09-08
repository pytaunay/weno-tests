function [qavep, qaven] = roe_average( q, q0, pbc )

%% Calculate q at i+1/2
rhol = q(:,1);
q2l = q(:,2);
q3l = q(:,3);

rhor = circshift(q(:,1),-1);
q2r = circshift(q(:,2),-1);
q3r = circshift(q(:,3),-1);


% Remove PBC if necessary
if(~pbc)
   rhor(end,1) = q0(end,1);   
   q2r(end,1) = q0(end,2);    
   q3r(end,1) = q0(end,3); 
    
end

ul = q2l./rhol;
e0l = q3l./rhol;

ur = q2r./rhor;
e0r = q3r./rhor;

sl = sqrt(rhol);
sr = sqrt(rhor);

% Average calculation
rhoave = sqrt(rhol.*rhor);
uave = (sl.*ul + sr.*ur)./(sl + sr);
e0ave = (sl.*e0l + sr.*e0r)./(sl+sr);

% Output at i+1/2
qavep = [rhoave,uave,e0ave];

%% Calculate q at i-1/2
rhol = circshift(q(:,1),1);
q2l = circshift(q(:,2),1);
q3l = circshift(q(:,3),1);

rhor = q(:,1);
q2r = q(:,2);
q3r = q(:,3);


% Remove PBC if necessary
if(~pbc)
   rhol(1,1) = q0(1,1);   
   q2l(1,1) = q0(1,2);    
   q3l(1,1) = q0(1,3); 
    
end

ul = q2l./rhol;
e0l = q3l./rhol;

ur = q2r./rhor;
e0r = q3r./rhor;

sl = sqrt(rhol);
sr = sqrt(rhor);

% Average calculation
rhoave = sqrt(rhol.*rhor);
uave = (sl.*ul + sr.*ur)./(sl + sr);
e0ave = (sl.*e0l + sr.*e0r)./(sl+sr);

% Output at i-1/2
qaven = [rhoave,uave.*rhoave,e0ave.*rhoave];


end

