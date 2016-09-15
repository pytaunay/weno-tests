function alpha = eulerMaxEigVal( qm,qp,GAM )

    am = speedOfSound(qm,GAM,'state');
    ap = speedOfSound(qp,GAM,'state');
    
    lam = zeros(size(qm,1),6);
    
    um = qm(:,2)./qm(:,1);
    
    up = qp(:,2)./qp(:,1);
    
    
    lam(:,1) = abs(um);
    lam(:,2) = abs(um - am);
    lam(:,3) = abs(um + am);
    
    lam(:,4) = abs(up);
    lam(:,5) = abs(up - ap);
    lam(:,6) = abs(up + ap);
 
    
    alpha = max(lam,[],2);

end

