%%% Function: uinit
%%% Returns the initial distribution of points to advect

function u0 = uinit( xcell )
  
    u0 = exp(-log(2)*(xcell+0.7).^2/9e-4).*(xcell>=-0.8).*(xcell<=-0.6);   
    u0 = u0 + 1.*(xcell >= -0.4).*(xcell <= -0.2);   
    u0 = u0 + (1-abs(10*xcell-1)).*(xcell>=0).*(xcell<=0.2);
    u0 = u0 + sqrt(1-100*(xcell-0.5).^2).*(xcell>=0.4).*(xcell<=0.6);
    

end

