function [Ts Rs Fr]=calculteFaraday(geometry,epsa,epsb,eps1,eps3,a1,a2,Rx,Ry...
  ,d1,d2,Na,nGx,nGy,k1,p,plotFT,plotWave,colorAng,theta,fi,rec)


nk=size(epsb,2);


d=d1+d2;

triang=0;

global MM;

global Kapa;
global bE;

w2c2=(k1)^2;

L=(Na-1)*a2+2*(Ry+d);

Lx1=2*nGx+1;
Lx2=3*nGx+2;

Ly1=2*nGy+1;


numG=(Lx1)*nGy;

dd=Lx1;

global inc_mode;
if(inc_mode==0)
E0=[0 0 1]';
else
E0=[1 0 0]';
end

if(nk==1) 
E0=[1];
end


kx=k1*sind(theta);
k1y=k1*cosd(theta);
k3=sqrt(eps3/eps1*w2c2);
bx=2*pi/a1;

by=pi/L;

pph=2/pi;

    dimx=nk*(numG+2*(Lx1));



if(p==1)
    disp('Computing Fourier series ..');
    Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);
    if(triang)
        FillKapaTriang(nGx,nGy,epsa,epsb,L,R,Na,a1,a2);  %triangular
        %lattivce not implemented properly.
    elseif (geometry==0 && Rx==Ry)
        if(nk==3)
         FillKapaCylinderAntiSymDef(nGx,nGy,epsa,epsb,L,Rx,Na,a1,a2,d1); 
       else
        FillKapaCylinderSym(nGx,nGy,epsa,epsb,L,Rx,Na,a1,a2,d1);
        end
    else
    if(geometry==1 && rec )
    FillKapaRectangleAntiSym(nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi);

  else

       FillKapaAntiSymNum(geometry,nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi);
 end
 end


    
    if(plotFT)
        
        
        ndx=30;
 	      ndy=40*Na;
        
        x=linspace(-a1/2,a1/2,ndx);
        y=linspace(-L,L,ndy);
        
        [xx yy]=ndgrid(x,y);
        
        zz=zeros(size(xx));
        
        for ix=1:ndx
            x1=x(ix);
            for iy=1:ndy

                y1=y(iy);
                tt=0;
                for n=-nGx:nGx
                    for m=-nGy:nGy
                        Gn=bx*n;
                        Gm=by*m;

                        TT=Kapa(n+2*nGx+1,m+2*nGy+1,4);
                        
                        tt=tt+TT*exp(1i*(Gn*x1+Gm*y1));
                    end
                end
                
                zz(ix,iy)=real(tt);
            end
        end
        
        
    %    X=abs(Kapa(:,:,1));

        
        figure(4);

        surf(xx,yy,zz);

        %  surf(X);
        axis equal;
        set(gca,'DataAspectRatio',[1 1 .05]);
%        return;
    end


    disp('Computing matrix, step 1...');
    
    
   
    MM=  zeros(dimx,dimx);
    

        
    bE=zeros(dimx,1);
    
    countG=0;
    
     
    nN=nGx;
    nM=nGy;
    for n=-nN:nN
        for m=1:nM
            
            countG=countG+1;
            countG1=0;
            
            for n1=-nN:nN
                
                dnn1=Lx1+n-n1;
                n1p=nN+1+n1;
                n1pp=Lx2+n1;
   
                kxn1=kx+n1*bx;
               % kxn1=0;
                
                kxn2=kxn1^2;

                
                for m1=1:nM
                    kym=m1*by;
                    dmm1=Ly1+(m-m1);
                    smm1=Ly1+m+m1;
                    countG1=countG1+1;    
                    
                    if(nk==1)
                     
                    Kapas(1:nk)=Kapa(dnn1,smm1,1:nk);                                 
                    Kapad(1:nk)=Kapa(dnn1,dmm1,1:nk);
                    
                      KapTenss=[Kapas(1) 0 0; 0 Kapas(1) 0 ;
                        0 0  Kapas(1)];

                      KapTensd=[Kapad(1) 0 0; 0 Kapad(1) 0 ;
                        0 0  Kapad(1)];
                    else
                    
                    Kapas(1:4)=Kapa(dnn1,smm1,1:4);                                 
                    Kapad(1:4)=Kapa(dnn1,dmm1,1:4);
                    
                    KapTenss=[Kapas(1) 0 1i*Kapas(4); 0 Kapas(2) 0 ;
                        -1i*Kapas(4) 0  Kapas(3)];

                    
                    KapTensd=[Kapad(1) 0 1i*Kapad(4); 0 Kapad(2) 0 ;
                        -1i*Kapad(4) 0  Kapad(3)];
                    end

                     Wp=zeros(3,3);
                    
                    
                    Wp(1,1)=(m1*by)^2;
                    Wp(1,2)=-n1*bx*m1*by;
                    Wp(2,1)= Wp(1,2);
                    Wp(2,2)=(n1*bx)^2;
                    
                    Wp(3,3)=Wp(1,1)+ Wp(2,2);
                    
                    Wm=zeros(3,3);
                    
                    Wm(1,1)= Wp(1,1);
                    Wm(1,2)=-Wp(1,2);
                    Wm(2,1)=-Wp(1,2);
                    Wm(2,2)=Wp(2,2);
                    Wm(3,3)=Wp(3,3);


                    FF=-(KapTensd*Wm-KapTenss*Wp);
                    

           
                    fRm1=pph/m1;  
                       
                      if(mod(m1,2)==0) 
                        fTm1=pph/m1;                    
                      else
                      fTm1=-pph/m1;
                    end
                    
                    
                     WT=zeros(3,3);
                     
                    WT(1,2)=1i*kxn1;
                    WT(2,1)= WT(1,2);
                    WT(2,2)=kxn2;
                      WT(3,3)=WT(2,2);
                   FFT=fTm1*(KapTensd+KapTenss)*WT;

                    
                    WR=zeros(3,3);

                    WR(1,2)=-1i*kxn1;
                    WR(2,1)= WR(1,2);
                    WR(2,2)=kxn2;
                   WR(3,3)=WR(2,2);

                    FFR=-fRm1*(KapTensd*WR+KapTenss)*WR;

                   
                    r1=(countG-1)*nk;
                    c1=(countG1-1)*nk;
                    for j=1:nk
                        for k=1:nk
                        
                        MM(r1+j,c1+k)= FF(j,k);
                                       
                         %% for Tn  
                         
                         i1=(countG-1)*nk+j;
                         i2=(numG+n1p-1)*nk+k;
                        MM((countG-1)*nk+j,(numG+n1p-1)*nk+k)= MM((countG-1)*nk+j,(numG+n1p-1)*nk+k)+FFT(j,k);
                         %% for Rn  
                        MM((countG-1)*nk+j,(numG+dd+n1p-1)*nk+k)= MM((countG-1)*nk+j,(numG+dd+n1p-1)*nk+k)+FFR(j,k);

                        end
                    end
                    
                    
                end
                
                
                               
                if(n1==0)
                    v=-FFR*E0;
                    for k=1:nk
                        bE((countG-1)*nk+k)=bE((countG-1)*nk+k)+v(k);
                    end
                    
                end
                   
                
                if(n==n1)
                    for k=1:nk
                        if(k~=2)
                            
                            MM((n1p+numG-1)*nk+k,((n1p-1)*nGy+m-1)*nk+k)=pi*m;
                        end
                    end
                    
                    if(mod(m,2)==0)
                        for k=1:nk
                            if(k~=2)
                                MM((n1pp+numG-1)*nk+k,((n1p-1)*nGy+m-1)*nk+k)=pi*m;
                            end
                        end
                    else
                        for k=1:nk
                            if(k~=2)
                                MM((n1pp+numG-1)*nk+k,((n1p-1)*nGy+m-1)*nk+k)=-pi*m;
                            end
                        end
                    end
                    
                end
                
            end
            
            
        end
        
        
    end

end


bN=zeros(dimx,1);

disp('Computing matrix, step 2...');

    NN=zeros(dimx,dimx);

countG=0;

for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    Gxpp=Lx2+Gx;
    
   
    for Gy=1:nGy
        
        countG=countG+1;
        for k=1:nk
                NN((countG-1)*nk+k,(countG-1)*nk+k)= w2c2;
        end
        
        if(mod(Gy,2)==0)
            for k=1:nk
                     NN((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=-pph/Gy*w2c2;
            end
        else
            for k=1:nk
                    NN((countG-1)*nk+k,(numG+Gxp-1)*nk+k)=pph/Gy*w2c2;

            end
        end
        
        for k=1:nk

                NN((countG-1)*nk+k,(numG+Gxpp-1)*nk+k)=pph/Gy*w2c2;
      
            if(Gx==0)
                
                bN((countG-1)*nk+k)=-pph/Gy*w2c2*E0(k);
                
            end
        end
        
    end
end


for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    Gxpp=Lx2+Gx;
    
    kxn1=kx+Gx*bx;
    kxn2=kxn1^2;
    
    
    if(k1>=abs(kxn1))
        krny=-sqrt(k1^2-kxn2);
    else
        krny=-1i*sqrt(kxn2-k1^2);
    end
    
    
    if(k3>=abs(kxn1))
        ktny=sqrt(k3^2-kxn2);
    else
        ktny=1i*sqrt(kxn2-k3^2);
    end
    

    
    for k=1:nk
        
       
            if(k~=2)
                NN((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
                
                
                NN((numG+Gxp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-(1i*L*krny+1);
                
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxp-1)*nk+k)=-(1i*L*ktny-1);
                
                NN((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=-1;
            else
                
                NN((numG+Gxp-1)*nk+k,(numG+Gxp-1)*nk+k)=1;
            
                NN((numG+Gxpp-1)*nk+k,(numG+Gxpp-1)*nk+k)=1;
            end
        
        if(Gx==0)
                      
            if(k~=2)
            bN((numG+Gxp-1)*nk+k)= (1i*L*k1y+1)*E0(k);
            bN((numG+Gxpp-1)*nk+k)=E0(k);
            end
        end
    end
end


bN=bN+bE;

  
disp('solving matrix...');


  NN=NN+MM;   

  x=linsolve(NN,bN);


Anm=zeros(2*nGx+1,nGy,nk);
Tn=zeros(2*nGx+1,nk);
Tn2=zeros(2*nGx+1,1);
Rn=zeros(2*nGx+1,nk);
Rn2=zeros(2*nGx+1,1);

ix=0;
for Gx=-nGx:nGx
    Gxp=nGx+1+Gx;
    for Gy=1:nGy
        for k=1:nk
            ix=ix+1;
            Anm(Gxp,Gy,k)=x(ix);
        end
    end
end


kp=numG*nk;
kpp=(numG+2*nGx+1)*nk;
for Gx=-nGx:nGx
    k=nGx+Gx;
   
    for j=1:nk
        Tn(k+1,j)=x(kp+k*nk+j);
        Tn2(k+1)= Tn2(k+1)+abs(Tn(k+1,j))^2;
        
        Rn(k+1,j)=x(kpp+k*nk+j);
        Rn2(k+1)= Rn2(k+1)+abs(Rn(k+1,j))^2;
    end
    

end


Ts=0;
Rs=0;

for Gx=-nGx:nGx
    k=nGx+1+Gx;
    kxn1=kx+Gx*bx;
    kxn2=kxn1^2;
    
    if(k1>=abs(kxn1))
        krny=-sqrt(k1^2-kxn2);
    else
        krny=-1i*sqrt(kxn2-k1^2);
    end
    
    
    if(k3>=abs(kxn1))
        ktny=sqrt(k3^2-kxn2);
    else
        ktny=1i*sqrt(kxn2-k3^2);
    end
        
    Ts=Ts+ktny/k3*sqrt(eps3/eps1)*Tn2(k)/cosd(theta);
        %Ts=Ts+ktny/k3*sqrt(eps3/eps1)*abs(Tn(k,3))^2/cos(thetad);

    
    Rs=Rs+abs(krny)/k1*Rn2(k)/cosd(theta);
end

Fr=0;


nL=max(int32(L*k1/2/pi)*17*2,35);



yy=linspace(0,L,nL);

Nx=1;

if(Nx==1)
    xx=zeros(1,1);
else
    xx=linspace(-a1/2,a1/2,Nx);
end

ff=zeros(Nx,nL,nk);
si=zeros(Nx,nL,nk);

for ix=1:Nx
    
    x=xx(ix);
    
    for k=1:nL
        y=yy(k);
        
        for Gx=-nGx:nGx
            
            kxn=kx+Gx*bx;
            
            Gxp=nGx+1+Gx;
            del=(Gx==0);
            for Gy=1:nGy
                for j=1:nk

		 ky0=Gy*by;

                    si(ix,k,j)=si(ix,k,j)+Anm(Gxp,Gy,j)*sin(ky0*y)*exp(1i*kxn*x);
                end
            end
            for j=1:nk
                ff(ix,k,j)=ff(ix,k,j)+(y*Tn(Gxp,j)+(L-y)*(Rn(Gxp,j)+del*E0(j)));
            end
        end
        
    end
end

E1=ff/L+si;

E2=real(E1);
E2(:,:,1);

%writeMeshAndField(Nx,nL,1,E2,2,Na);


if(plotWave)
    figure(5)
    
    set(gca,'DataAspectRatio',[1 1 1]);
    axis([-2 2 0 L -2 2]);
    az = 40;
    el = 30;
    %  az=90;
    %  el=0;
    view(az, el);
    
    
    
    hold on
    
    for ix=1:Nx
        x=xx(ix);
        Vx=E2(ix,1,1);
        Vy=E2(ix,1,2);
        Vz=E2(ix,1,3);
        
        
        for k=1:nL
            
            y=yy(k);
            Vxp=Vx;
            Vx=E2(ix,k,1);
            
            Vyp=Vy;
            Vy=E2(ix,k,2);
            
            Vzp=Vz;
            Vz=E2(ix,k,3);
            
            
            color = 'r';
            
            
            if(k>1)
                h1(k)=plot3([Vxp+x Vx+x],[Vyp+yy(k-1) Vy+y],[Vzp Vz],'LineWidth',2,'Color','b');
                h2(k)= plot3([-Vxp+x -Vx+x],[yy(k-1) y],[-Vzp -Vz],'LineWidth',2,'Color','b');
                h3= plot3([0 0],[0 L],[0 0],'LineWidth',2,'Color','k');
                
                
            end
            
            
            arrow1(k) = arrow3d([x x+Vx],[y y+Vy],[0 Vz],.92,.01,.02,color);
            arrow1(k) = arrow3d([x x-Vx],[y y-Vy],[0 -Vz],.92,.01,.02,color);
            
            
            
        end
        
    end;
end


if(nk==3)

if(inc_mode==0)
tn0=E2(1,1,1)/E2(1,1,3);
else
tn0=-E2(1,1,3)/E2(1,1,1);
end

ang0=atan(tn0);
nan=0;

tans1=zeros(nL,1);
for k=1:nL
    
    Vx=E2(1,k,1);
    Vz=E2(1,k,3);
    
    magn=sqrt(Vx^2+Vz^2);
    
    if( magn>-.1)
    if(inc_mode==0)
    tans1(k)=Vx/Vz;
else
    tans1(k)=-Vz/Vx;
end
        nan=nan+1;
    else
        tans1(k)=1e10;
    end
end

tans=zeros(nan,1);
y=zeros(nan,1);
j=0;
for k=1:nL
    if(tans1(k)~=1e10)
        j=j+1;
        
        tans(j)=tans1(k);
        
        y(j)=yy(k);
    end
end;

angs=zeros(nan,1);
angs(1)=atan(tans(1));

for k=2:nan
    
    tnp=tans(k-1);
    tn=tans(k);
    angs(k)=-atan(tn)*180/pi;

    if(angs(k)>0)
  
%    angs(k)=-180+angs(k);
    
    end

    %  angs(k)=E2(1,k,3);
    % angs(k)=angs(k-1)+atan((tn-tnp)/(1+tn*tnp))*180/pi;
    
    
end



Fr=(angs(nan)-ang0);

end
%Eout=sqrt(E2(1,nL,1)^2+E2(1,nL,2)^2+E2(1,nL,3)^2);

if(plotWave)
    
    figure(7)
    
    plot(y,angs,colorAng);
    hold on
    
    %angsHomog=-0.5*y*k1*imag(epsb(1,3))/sqrt(epsb(1,1))*180/pi;
    %plot(y,angsHomog,'-k');
    
end

end

