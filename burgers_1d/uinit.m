%%% Function: uinit
%%% Returns the initial distribution of points to advect

function u0 = uinit( xcell )
  

%u0 = 0.25 + 0.5*sin(pi*xcell);

u0 = (xcell>=0) .* (xcell <= 0.5);

end

