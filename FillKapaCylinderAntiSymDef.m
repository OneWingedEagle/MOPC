function FillKapaCylinderAntiSymDef(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

disp('anti-sym Defect');
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

ff=pi*R*R/(a1*a2);


KapaUnit=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);


for dGx=-2*nGx:2*nGx
    Gn=2*pi*dGx/a1;
    for dGy=-2*nGy:2*nGy
        Gm=dGy*pi/L;
        GnmR=sqrt(Gn*Gn+Gm*Gm)*R;
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;
        
        if(dGx==0&& dGy==0)
            for k=1:nk
                  KapaUnit(dGxp,dGyp,k)=(invepsb(k)+ff*(invepsa(k)-invepsb(k)))/(2*Na);
                  
                  Kapa0Defect(1,1,k)=invepsb(k)/(2*Na);
            end
            
         else
            tt=ff*besselj(1,GnmR)/GnmR;
            
            for k=1:nk

                    KapaUnit(dGxp,dGyp,k)=tt*(invepsa(k)-invepsb(k))/(Na);;
            end
            
            
        end
        
        
    end
end

global ndef


by=pi/L;


KapaDefect=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);
 


if(ndef>0)
 
     for dGy=-2*nGy:2*nGy
        
        dGyp=dGy+1+2*nGy;

        tt=dGy*by*a2/2;
          
          for k=1:nk
                if(dGy==0)
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L);
                else
                  KapaDefect(1,dGyp,k)=invepsb(k)*a2/(2*L)*sin(tt)/(tt);
             end
         
         end
    
        
    end
 
end
 
Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);



by=pi/L;
for n=-Na:Na-1 


 for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle=exp(-1i*(n+.5)*by*dGy*a2);
        
   if(n>=-ndef && n<ndef)

    for k=1:nk
       if(k!=4 || n>=0)
       Kapa(1,dGyp,k)=  Kapa(1,dGyp,k)+KapaDefect(1,dGyp,k)*twindle;
        else
       Kapa(1,dGyp,k)=  Kapa(1,dGyp,k)-KapaDefect(1,dGyp,k)*twindle;
        end
     end
   else
     for k=1:nk
       if(k!=4 || n>=0)
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)+KapaUnit(:,dGyp,k)*twindle;
        else
         Kapa(:,dGyp,k)=Kapa(:,dGyp,k)-KapaUnit(:,dGyp,k)*twindle;
        end
     end
   
   end
        
        
    end
  end

  
end



