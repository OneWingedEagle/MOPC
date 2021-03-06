function FillKapaAntiSymNumOld(geometry,nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi)

nk=4;
if(length(epsa)==1)
nk=1;
end

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

end

global Kapa;

 Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);


ngridx=40;
##if(nGx==0) 
##ngridx=30;
##end
ngridy=40;
Lh=Na*a2-a2/2;

invep=zeros(ngridx,Na*2*ngridy,nk);
rotMat=[cosd(fi) sind(fi);-sind(fi) cosd(fi)];


xyr=[0  0]';
%The following loop carries out the definition of the unit cell.
for n=0:2*Na-1
    nx=1;
    for countX=-a1/2:a1/ngridx:a1/2
        ny=1;
        for countY=-a2/2:a2/ngridy:a2/2
            %Following condition allows to define the circle with of radius r
			xy=[countX countY]';
			xyr=rotMat*xy;
			countXr=xyr(1);
           	countYr=xyr(2);
				
            if(geometry==0)
                inside=(countXr/Rx)^2+(countYr/Ry)^2<1;
            else
                 inside =false;
                  if(Rx>0 && Ry>0) 
                     inside=abs(countXr/Rx)<=1 && abs(countYr/Ry)<=1;
                  end
            end
            
            if(inside)
                
                ip=invepsa;

            else
                
                ip=invepsb;
 
            end
            
            
            
            for k=1:nk
                inveps(nx,ny+ngridy*n,k)=ip(k);
            end
            
            xSet(nx)=countX;
            ySet(ny+ngridy*n)=countY+a2*n-Lh;
            ny=ny+1;
        end
        nx=nx+1;
    end
end

for n=1:length(ySet)/2
    inveps(:,n,nk)=-inveps(:,n,nk);
end

MtNt=(length(xSet)-1)*(length(ySet)-1);
%The next loop computes the Fourier expansion coefficients
bx=2*pi/a1;

by=pi/L;


for dGx=-2*nGx:2*nGx
    
    for dGy=-2*nGy:2*nGy
        
        
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        for nx=1:length(xSet)-1
            for ny=1:length(ySet)-1
                x=xSet(nx);
                y=ySet(ny);
                tt=dGx*bx*x+dGy*by*y;
                
                for k=1:nk
                    Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+inveps(nx,ny,k)*exp(-1i*(tt));
                end
            end
        end
        
        
    end
end

Kapa=  Kapa/MtNt;

end