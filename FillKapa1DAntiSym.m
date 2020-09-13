function FillKapa1DAntiSym(nGx,nGy,epsa,epsb,L,R,Na,a1,a2,d1)
global Kapa;

nk=4;
if(length(epsa)==1)
nk=1;
endif


invepsa1=inv(epsa);
invepsb1=inv(epsb);

if(nk==1)

invepsa=[invepsa1(1,1)];
invepsb=[invepsb1(1,1)];

else
invepsa=[invepsa1(1,1)  invepsa1(2,2) invepsa1(3,3) imag(invepsa1(1,3))];
invepsb=[invepsb1(1,1) invepsb1(2,2) invepsb1(3,3)  imag(invepsb1(1,3))];

endif

 w0=pi/L;

  for n=1:Na
  
    for dGy=-2*nGy:2*nGy
       
        Gm=dGy*pi/L;
        GnmR=sqrt(Gm*Gm)*R;
        dGyp=dGy+1+2*nGy;
        
        if(dGy==0)
         tt=R/L;
            for k=1:nk
                if(k~=4)
                    Kapa(1,dGyp,k)= Kapa(1,dGyp,k)+exp(-1i*w0*n*a2)*tt*(invepsa(k)-invepsb(k));
                 else
                     Kapa(1,dGyp,k)=Kapa(1,dGyp,k)+exp(-1i*w0*n*a2)*tt*tt*(invepsa(k)-invepsb(k));
                endif
            end         
        else
            
          tt=1./(2*L)*sin(dGy*pi*R/(2*L))/(dGy*pi/(2*L));
            
            for k=1:nk
                if(k~=4)
                    Kapa(1,dGyp,k)=Kapa(1,dGyp,k)+exp(-1i*w0*n*a2)*tt*tt*(invepsa(k)-invepsb(k));
                 else
                    Kapa(1,dGyp,k)=Kapa(1,dGyp,k)+exp(-1i*w0*n*a2)*tt*tt*(invepsa(k)-invepsb(k));
               endif
            end
            
        endif
                       
        end
      end
          
end

