% ORIGINAL SOURCE: https://www3.nd.edu/~gtryggva/
% https://www3.nd.edu/~gtryggva/CFD-Course2017/Project-3-2017.pdf
% https://www3.nd.edu/~gtryggva/CFD-Course2017/Euler2D-ZB-upw-Oblique.m

nx=256; ny=64; xl=4.0;yl=1; tfinal=0.2; time=0; gg=1.4;
p_left=116.5;p_right=1;
r_left=8;r_right=1.4;
u_left=8.25*cos(pi/6); u_right=0.0;
v_left=-8.25*sin(pi/6); v_right=0.0;

r=zeros(nx,ny);p=zeros(nx,ny);rE=zeros(nx,ny); c=zeros(nx,ny);
ru=zeros(nx,ny);u=zeros(nx,ny);mx=zeros(nx,ny);
rv=zeros(nx,ny);v=zeros(nx,ny);my=zeros(nx,ny);
F1=zeros(nx,ny);F2=zeros(nx,ny);F3=zeros(nx,ny);F4=zeros(nx,ny);
G1=zeros(nx,ny);G2=zeros(nx,ny);G3=zeros(nx,ny);G4=zeros(nx,ny);
h=xl/(nx-1);for i=1:nx,for j=1:ny,x(i,j)=h*(i-1);y(i,j)=h*(j-1);end,end;

for i=1:nx,for j=1:ny,r(i,j)=r_right;
    r(i,j)=r_right;ru(i,j)=0.0; rv(i,j)=0; rE(i,j)=p_right/(gg-1); end;end

xo=1/6;    
for i=1:nx,for j=1:ny,
  if x(i,j) <  xo+(y(i,j)/sqrt(3));
   r(i,j)=r_left; rE(i,j)=p_left/(gg-1)+0.5*r_left*(u_left^2 + v_left^2); 
   ru(i,j)=r_left*u_left; rv(i,j)=r_left*v_left;
  end
end;end

cmax=sqrt( max(gg*p_right/r_right,gg*p_left/r_left) ); 
dt=0.145*h/cmax; maxstep=tfinal/dt;

for istep=1:maxstep
  for i=1:nx,for j=1:ny, p(i,j)=(gg-1)*(rE(i,j)-...
               0.5*((ru(i,j)*ru(i,j)+rv(i,j)*rv(i,j))/r(i,j)));end,end
  for i=1:nx,for j=1:ny,c(i,j)=sqrt( gg*p(i,j)/r(i,j) );end,end
  for i=1:nx,for j=1:ny,u(i,j)=ru(i,j)/r(i,j);v(i,j)=rv(i,j)/r(i,j);end,end;
  for i=1:nx,for j=1:ny,mx(i,j)=u(i,j)/c(i,j);my(i,j)=v(i,j)/c(i,j);end,end;

  for i=1:nx-1, for j=1:ny-1					% Find fluxes
    F1(i,j)=0.5*(ru(i+1,j)+ru(i,j))-0.5*(abs(ru(i+1,j))-abs(ru(i,j)));

    F2(i,j)=0.5*(u(i+1,j)*ru(i+1,j)+p(i+1,j)+u(i,j)*ru(i,j)+p(i,j))...
         -0.5*(abs(u(i+1,j))*ru(i+1,j)-abs(u(i,j))*ru(i,j))...
         -0.5*(p(i+1,j)*mx(i+1,j)-p(i,j)*mx(i,j));

    F3(i,j)=0.5*(u(i+1,j)*rv(i+1,j)+u(i,j)*rv(i,j))...
         -0.5*(abs(u(i+1,j))*rv(i+1,j)-abs(u(i,j))*rv(i,j));
         
    F4(i,j)=0.5*(u(i+1,j)*(rE(i+1,j)+p(i+1,j))+u(i,j)*(rE(i,j)+p(i,j)))...
         -0.5*(abs(u(i+1,j))*rE(i+1,j)-abs(u(i,j))*rE(i,j))...
         -0.5*(p(i+1,j)*c(i+1,j)-p(i,j)*c(i,j));

    if mx(i,j) > 1,  F2(i,j)=ru(i,j)*u(i,j)+p(i,j);       
                           F4(i,j)=(rE(i,j)+p(i,j))*u(i,j);end
    if mx(i,j) < -1, F2(i,j)=ru(i+1,j)*u(i+1,j)+p(i+1,j); 
                      F4(i,j)=(rE(i+1,j)+p(i+1,j))*u(i+1,j);end

    G1(i,j)=0.5*(rv(i,j+1)+rv(i,j))-0.5*(abs(rv(i,j+1))-abs(rv(i,j)));

    G2(i,j)=0.5*(v(i,j+1)*ru(i,j+1)+v(i,j)*ru(i,j))...
         -0.5*(abs(v(i,j+1))*ru(i,j+1)-abs(v(i,j))*ru(i,j));

    G3(i,j)=0.5*(v(i,j+1)*rv(i,j+1)+p(i,j+1)+v(i,j)*rv(i,j)+p(i,j))...
         -0.5*(abs(v(i,j+1))*rv(i,j+1)-abs(v(i,j))*rv(i,j))... 
         -0.5*(p(i,j+1)*my(i,j+1)-p(i,j)*my(i,j));

    G4(i,j)=0.5*(v(i,j+1)*(rE(i,j+1)+p(i,j+1))+v(i,j)*(rE(i,j)+p(i,j)))...
         -0.5*(abs(v(i,j+1))*rE(i,j+1)-abs(v(i,j))*rE(i,j))...
         -0.5*(p(i,j+1)*c(i,j+1)-p(i,j)*c(i,j));

    if my(i,j) > 1,  G3(i,j)=rv(i,j)*v(i,j)+p(i,j);       
                           G4(i,j)=(rE(i,j)+p(i,j))*v(i,j);end
    if my(i,j) < -1, G3(i,j)=rv(i,j+1)*v(i,j+1)+p(i,j+1); 
                      G4(i,j)=(rE(i,j+1)+p(i,j+1))*v(i,j+1);end                
  end, end


  for i=2:nx-1 ; for j=2:ny-1   				% Update solution
    r(i,j) =r(i,j) -(dt/h)*(F1(i,j)-F1(i-1,j)+G1(i,j)-G1(i,j-1));
    ru(i,j)=ru(i,j)-(dt/h)*(F2(i,j)-F2(i-1,j)+G2(i,j)-G2(i,j-1));
    rv(i,j)=rv(i,j)-(dt/h)*(F3(i,j)-F3(i-1,j)+G3(i,j)-G3(i,j-1));
    rE(i,j)=rE(i,j)-(dt/h)*(F4(i,j)-F4(i-1,j)+G4(i,j)-G4(i,j-1));
  end; end
  
  % Bottom
  no=floor(nx*xo/xl)+2;
  r(no:nx-2,1)=r(no:nx-2,2);ru(no:nx-2,1)=ru(no:nx-2,2);rE(no:nx-2,1)=rE(no:nx-2,2);
  
  %Top
  xs=xo+((1+20*time))/sqrt(3); ns=floor(nx*xs/xl)+1;
  for i=2:ns,
    r(i,ny)=r_left;  ru(i,ny)=r_left*u_left; rv(i,ny)=r_left*v_left;
    rE(i,ny)=p_left/(gg-1)+0.5*(ru(i,ny)*ru(i,ny)/r(i,ny) + rv(i,ny)*rv(i,ny)/r(i,ny));
  end
  
%  r(1,2:ny-2)=r(2,2:ny-2);rv(1,2:ny-2)=rv(2,2:ny-2);rE(1,2:ny-2)=rE(2,2:ny-2);
  r(nx,2:ny-2)=r(nx-1,2:ny-2);rv(nx,2:ny-2)=rv(nx-1,2:ny-2);rE(nx,2:ny-2)=rE(nx-1,2:ny-2);

  time=time+dt, istep
  Pplot = 
  contour(x,y,r,40); axis equal, axis([0, xl, 0, yl]);pause(0.001)
end
% plot(x,r,'k','linewidth',2);hold on
% set(gca,'Box','on'); set(gca,'Fontsize',24, 'LineWidth',2)
% text(5,0.9,'Density','Fontsize',24)
% print -depsc ZBResults1     mesh(x,y,r)    quiver(x,y,u,v)

axis([0, xl*3/4, 0, yl])
