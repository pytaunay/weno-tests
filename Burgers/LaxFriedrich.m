function F = LaxFriedrich(UL,UR)

% Input: 
% - QL: value at the left side of the boundary 
% - Qr: value at the right side of the boundary 
% [       QL][QR     ]...

% alpha = max f'(u), where f = u^2/2 for Burgers

alpha = max(abs(UL),abs(UR));

F = 1/2*(BurgersFlux(UL)+BurgersFlux(UR)-alpha*(UR-UL));

end

