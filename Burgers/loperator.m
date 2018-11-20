function URET = loperator(U,DX,UL,UR,NCELL,flux)
% Calculates the "L" operator in the linearized hyperbolic equation:
% du/dt = -1/dx*(F(i+1/2)-F(i-1/2))

% Requires reconstruction of Q at cell boundaries to calculate the flux at
% the boundaries
URET = zeros(NCELL,1);

U1 = 0;
U2 = 0;
U3 = 0;
U4 = 0;

for i = 1:NCELL

   if( i == 1)

       U1 = UL;            % 0, Right
       U2 = UL;            % 1, Left
       U3 = WENO5(U,i,'L',NCELL);
       U4 = WENO5(U,i+1,'R',NCELL);   
       
   elseif( i == NCELL)
       
       U1 = U3;     % i-1, right
       U2 = U4;     % i, left       
       U3 = UR;
       U4 = UR;
       
   else
       U1 = U3;     % i-1, right
       U2 = U4;     % i, left   
       U3 = WENO5(U,i,'L',NCELL);
       U4 = WENO5(U,i+1,'R',NCELL);             
   end

   if( strcmp(flux,'LF') )
       URET(i,1) = LaxFriedrich(U3,U4);
       URET(i,1) = URET(i,1) - LaxFriedrich(U1,U2);
   elseif ( strcmp(flux,'GOD') )
       URET(i,1) = Godunov(U3,U4);
       URET(i,1) = URET(i,1) - Godunov(U1,U2);
   end
   
   URET(i,1) = -URET(i,1)/DX;
     
end

end