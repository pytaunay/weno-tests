function alpha = eulerMaxEigVal( qm,qp,GAM )

    am = speedOfSound(qm,GAM);
    ap = speedOfSound(qp,GAM);
    
    lam = zeros(size(qm,1),6);
    

    um = qm(:,2)./qm(:,1);
    up = qp(:,2)./qp(:,1);
    
    
    lam(:,1) = um;
    lam(:,2) = um -am;
    lam(:,3) = um + am;
    
    lam(:,4) = up;
    lam(:,5) = up - ap;
    lam(:,6) = up + ap;
    
    alpha = max(lam,[],2);


end

