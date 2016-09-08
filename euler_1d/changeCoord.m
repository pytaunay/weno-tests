function r = changeCoord( R,q )

r = zeros(size(q,1),3);

for i = 1:size(q,1)
    qi = q(i,:)';
    Ri = R{i};
    
    r(i,:) = (Ri*qi)';

end

