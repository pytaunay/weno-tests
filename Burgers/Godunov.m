function F = Godunov(UL,UR)

% Input: 
% - QL: value at the left side of the boundary 
% - Qr: value at the right side of the boundary 
% [       QL][QR     ]...

% f = u^2/2 for Burgers
if(UL <= UR)
    F = min(UL^2/2,UR^2/2);
else
    F = max(UL^2/2,UR^2/2);
end


end

