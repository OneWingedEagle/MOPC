function FillKapaCylinderAntiSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

nk=length(epsa);

invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

endif

ff=Na*pi*R*R/(a1*L);
d11=d1/L;

for dGx=-2*nGx:2*nGx
    Gn=2*pi*dGx/a1;
    for dGy=-2*nGy:2*nGy
        Gm=dGy*pi/L;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                if(k~=4)
                    Kapa(dGxp,dGyp,k)=1*(ff*invepsa(k)+(1-ff)*invepsb(k))+2*d11*(invepsa(k)-invepsb(k));
                end
                
            end
            
        elseif(dGy==0)
            tt=2*ff*besselj(1,GnmR)/GnmR;
            
            for k=1:nk
                if(k~=4)
                    Kapa(dGxp,dGyp,k)=tt*(invepsa(k)-invepsb(k));
                end
            end
            
        else
            tt=a2*dGy*pi/(2*L);
            
            
            for k=1:nk
                if(k==4)
                    factm=-1i*sin(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*factm*vv*(invepsa(k)-invepsb(k));
                    if(dGx==0)
                        vvL=2*sin(dGy*pi*.5*(L-2*d1)/L)/(dGy*pi);
                        factm=-1i*sin(dGy*pi*.5);
                        
                        Kapa(dGxp,dGyp,k)=Kapa(dGxp,dGyp,k)+factm*vvL*invepsb(k);
                    end
                else
                    factm=cos(dGy*pi/2)*sin(Na*tt)/sin(tt)/Na;
                    vv=2*ff*factm*besselj(1,GnmR)/GnmR;
                    
                    Kapa(dGxp,dGyp,k)=1*vv*(invepsa(k)-invepsb(k));
                end
                
                
                if(dGx==0 && d11>0 && k~=4)
                    aa= pi*dGy*d11;
                    
                    zz=1*sin(aa)/aa*d11+2*cos(pi*dGy*(L-d1/2)/L)*sin(aa/2)/aa*d11;
                    Kapa(dGxp,dGyp,k)= Kapa(dGxp,dGyp,k)+zz*(invepsa(k)-invepsb(k));
                end
                
            end
            
        end
        
        
    end
end
end


