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

w0x=2*pi/a1;
w0y=pi/L;

for dGx=-2*nGx:2*nGx
   
   for dGy=-2*nGy:2*nGy
        
        
        dGxp=dGx+1+2*nGx;
        dGyp=dGy+1+2*nGy;

        four_coef=0;
          if(dGx==0 && dGy==0)
            four_coef=rect_size/cell_size;
         elseif(dGx==0)
           four_coef=fillx*sin(dGy*w0y*Ry)/(dGy*w0y);
         elseif(dGy==0)
            four_coef=filly*sin(dGx*w0x*Rx)/(dGx*w0x);
       else
            four_coef=sin(dGx*w0x*Rx)*sin(dGy*w0y*Ry)/(dGx*dGy*w0x*w0y);
        end
        
 %         four_coef=four_coef*exp(-1i*w0y*a2);
        for n=-(Na):Na-1
       %   four_coef=four_coef+four_coef*exp(-1i*n*w0y*a2);
          end
        
           for k=1:nk
             
             if(dGx==0 && dGy==0)
                Kapa(dGxp,dGyp,k)==four_coef*invepsa(k)+(1-four_coef)*invepsb(k);
           else
               Kapa(dGxp,dGyp,k)=four_coef*(invepsa(k)-invepsb(k));

           end
           
               
               
           end
        

  end

end

end

