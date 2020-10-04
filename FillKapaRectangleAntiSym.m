function FillKapaRectangleAntiSym(nGx,nGy,epsa,epsb,L,Rx,Ry,Na,a1,a2,d1,fi)
global Kapa;

disp("rectangle anti-sym");
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

cell_size=a1*a2*Na*2;
rect_size=2*Rx*2*Ry;

fillx=2*Rx/a1;
filly=Ry/L;

bx=2*pi/a1;
by=pi/L;


KapaUnit=zeros(1,4*nGy+1,4)+1i*zeros(1,4*nGy+1,4);



for dGx=-2*nGx:2*nGx
      dGxp=dGx+1+2*nGx;
   ttx=dGx*bx*Rx;
   for dGy=-2*nGy:2*nGy
        
        tty=dGy*by*Ry;
     
        dGyp=dGy+1+2*nGy;
  for k=1:nk
      if(dGx==0)
       four_coefx= (invepsb(k)+fillx*(invepsa(k)-invepsb(k)));
        else
          four_coefx=(invepsa(k)-invepsb(k))*2*Rx/(a1)*sin(ttx)/(ttx);
       end
        
        
         if(dGy==0)
       four_coefy= (invepsb(k)+fillx*(invepsa(k)-invepsb(k)))/(2*Na);
        else
          four_coefy=(invepsa(k)-invepsb(k))*2*Ry/(a1)*sin(tty)/(tty)/(2*Na);
       end
       
        KapaUnit(dGxp,dGyp,k)=four_coefx*four_coefy;


  end

end

end

Kapa=zeros(4*nGx+1,4*nGy+1,4)+1i*zeros(4*nGx+1,4*nGy+1,4);


global ndef;
global defstart;

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

ndef
defstart

defcount=0;
isDef=zeros(2*Na,1);
for n=0:Na-1 
  np=n+Na+1;
 if(n+1>=defstart && defcount<ndef)
 defcount=defcount+1;
 isDef(np,1)=1;
 isDef(Na-n,1)=1;
 endif
end


by=pi/L;
for n=-Na:Na-1 
 np=n+Na+1;

 for dGy=-2*nGy:2*nGy
       
      dGyp=dGy+1+2*nGy;
        twindle=exp(-1i*(n+.5)*by*dGy*a2);
        
     if(isDef(np,1))

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
Kapa;
  
end




