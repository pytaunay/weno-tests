function r = changeCoord( R,q )

r = zeros(size(q,1),3);

for i = 1:size(q,1)
    % Change from row to column representation
    qi = q(i,:)';
    Ri = R{i};
    
    % Go back to row representation
    r(i,:) = (Ri*qi)';

end

