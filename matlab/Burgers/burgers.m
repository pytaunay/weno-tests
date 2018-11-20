% Solve the inviscid Burgers equation and compares it to the analytical solution
% FVM - based solution with:
% TVD RK-3 (Shu)
% Lax-Friedrich Flux 
% WENO order 5 for flux reconstruction

clear all

% CONSTANTS AND INITIAL VALUES
NRUN = 4;

L1Plot = zeros(NRUN,1);
L2Plot = zeros(NRUN,1);
LinfPlot = zeros(NRUN,1);
NCELLPlot = zeros(NRUN,1);
j = 1;

figure()
for NCELL = 100:100:100;
    
    NCELLPlot(j,1) = NCELL;
    
    XMAX = 1;
    XMIN = -1;
    DX = (XMAX-XMIN)/NCELL;
    CFL = 0.4;
    
    TMAX = 2;
    
    U = zeros(NCELL,1);
    XCELL = zeros(NCELL,1);
    
    
    UL = 1;
    UR = -0.5;
        
    s = 1/2*(UR+UL);
    
    %%% GODUNOV
    flux = 'GOD';
    
    for i = 1:NCELL
        
        % CELL CENTER
        XCELL(i,1) = (i-1)*DX+DX/2+XMIN;
        if( XCELL(i,1) < 0 )
            U(i,1) = UL;
        else
            U(i,1) = UR;
        end
    end
    
    T = 0;
    Uk = zeros(NCELL,1);
    TV = zeros(int64(TMAX/0.001),1);
    l = 1;
    while (T<TMAX)
                
        % Calculate DT
%         DT = CFL*DX*max(abs(U(:,1))+s);
        DT = CFL*DX/max(abs(U(:,1)));
        UP = U;
        
        % Check we indeed have TVD
        for i = 2:NCELL
           TV(l,1) = TV(l,1) + abs(UP(i,1)-UP(i-1,1)); 
        end
        l = l+1;
        
        % TVD RK3
        Uk  = U         +   DT*loperator(U,DX,UL,UR,NCELL,flux);
        Uk  = 0.75*U    +   0.25*Uk +   0.25*DT*loperator(Uk,DX,UL,UR,NCELL,flux);
        U   = 1/3*U     +   2/3*Uk  +   2/3*DT*loperator(Uk,DX,UL,UR,NCELL,flux);
        
        % Simple Euler
%         U = U + DT*loperator(U,DX,UL,UR,NCELL,flux);
        
        
        
        T = T+DT;
        

    end
    
    
    
    subplot(NRUN,1,j);
    plot(XCELL(:,1),U(:,1),'--b.');
    UGOD = U;
    
    
%     %%% LAX FRIEDRICH
%     flux = 'LF';
%     T = 0;
%     Uk = zeros(NCELL,1);
%     for i = 1:NCELL
%         
%         % CELL CENTER
%         XCELL(i,1) = (i-1)*DX+DX/2+XMIN;
%         if( XCELL(i,1) < 0 )
%             U(i,1) = UL;
%         else
%             U(i,1) = UR;
%         end
%     end
%     
%     while (T<TMAX)
%         
%         
%         % Calculate DT
%         DT = CFL*DX*max(abs(U(:,1))+s);
%         
%         % TVD RK3
%         Uk  = U         +   DT*loperator(U,DX,UL,UR,NCELL,flux);
%         Uk  = 0.75*U    +   0.25*Uk +   0.25*DT*loperator(Uk,DX,UL,UR,NCELL,flux);
%         U   = 1/3*U     +   2/3*Uk  +   2/3*DT*loperator(Uk,DX,UL,UR,NCELL,flux);
%         
%         T = T+DT;
%         
%         % Check we indeed have TVD
%         
%         
%     end
%     
%     hold on
%     plot(XCELL(:,1),U(:,1),'--g.');
%     ULF = U;
    
    %%% EXACT SOLUTION
   
    XMINEX = XMIN+DX/2;
    XMAXEX = XMAX-DX/2;
%     DXEX = (XMAXEX-XMINEX)/NCELLEX;
    DXEX = DX/100;
    NCELLEX = (XMAXEX-XMINEX)/DXEX+1;
    NCELLEX = int64(NCELLEX);
    
    Xex = zeros(NCELLEX,1);
    Uex = zeros(NCELLEX,1);   
    
    s = 1/2*(UL+UR);
    
    for k = 1:NCELLEX
        
        %CELL CENTER
        Xex(k,1) = XMINEX+ double((k-1))*DXEX;

        if(UL > UR)
            if(Xex(k,1) < s*TMAX)
                Uex(k,1) = UL;
            else
                Uex(k,1) = UR;
            end
        elseif(UL < UR)
            if(Xex(k,1) < UL*TMAX)
                Uex(k,1) = UL;
            elseif( UL*TMAX < Xex(k,1) && Xex(k,1) < UR*TMAX )
                Uex(k,1) = Xex(k,1)/TMAX;
            else
                Uex(k,1) = UR;
            end
        end
        
    end
    
    hold on
    plot(Xex(:,1),Uex(:,1),'r');
%     legend('Godunov','Lax-Friedrich','Exact');
     legend('Godunov','Exact');

    % ERROR CALCULATION GODUNOV
    L1err = 0;
    L2err = 0;
    L2den = 0;
    Linferr = 0;
    Linfden = zeros(NCELL,1);
    abslist = zeros(NCELL,1);
    Ilist = zeros(NCELL,1);
    % % Error calculation
    for i = 1:NCELL
%         x = XCELL(i)-XMIN;
%         idx = x/DX;
        idx = (XCELL(i)-XMINEX)/DXEX+double(1);
        
        I = int64(idx);
        
        abslist(i,1) = abs(UGOD(i,1)-Uex(I,1));
%         L1err = L1err + abslist(i,1)/abs(Uex(I,1));
        L1err = L1err + abslist(i,1);
        L2err = L2err + abslist(i,1)^2;
        L2den = L2den + Uex(I,1)^2;
        Linfden = abs(Uex(I,1));
        Ilist(i,1) = I;
        
    end
    
    
%     L1Plot(j,1) = L1err/NCELL;
    L1Plot(j,1) = L1err;
    L2Plot(j,1) = sqrt(L2err/L2den);
    LinfPlot(j,1) = max(abslist./Linfden);
    j = j+1;

end


figure()
plot(NCELLPlot,L1Plot,'b');
hold on
plot(NCELLPlot,L2Plot,'r');
hold on
plot(NCELLPlot,LinfPlot,'g');
legend('L1','L2','Linf');

