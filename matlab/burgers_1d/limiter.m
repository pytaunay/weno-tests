function out = limiter( r, limiterType )

if(strcmp(limiterType,'superbee'))
    out = max(0, max( min(2*r,1), min(r,2)));
elseif( strcmp(limiterType,'vanLeer'))
    out = (r + abs(r))./(1+abs(r));
elseif( strcmp(limiterType,'koren'))
    out = max(0, min( 2*r,min( (2+r)/3 , 2)));
    
else
    out = max(0, min(r,1));
end

end

