function faraday

clear all

rotationGraphColor='-ob';
transmitGraphColor='-or';

[filename1,filepath1]=uigetfile('*.txt', 'Selectinput file')
 cd(filepath1)
 fid= fopen(filename1)
%fid = fopen('input.txt','rt');

line=getNewDataLine(fid);
numbs = str2num(line);
geometry=numbs(1);

line=getNewDataLine(fid);
numbs = str2num(line);
Rx=numbs(1);
Ry=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
fi=numbs(1);

line=getNewDataLine(fid);
numbs = str2num(line);
theta=numbs(1);

line=getNewDataLine(fid);
numbs = str2num(line);
ax=numbs(1);
ay=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
eps1=numbs(1);
eps3=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
ndata=length(numbs)

epsbx=numbs(1);
if(ndata>1)
epsby=numbs(2);
epsbz=numbs(3);
gamab=numbs(4);
end

line=getNewDataLine(fid);
numbs = str2num(line);
epsax=numbs(1);
if(ndata>1)
epsay=numbs(2);
epsaz=numbs(3);
gamaa=numbs(4);
end

if(ndata>1)

epsb=[epsbx 0 -1i*gamab;0 epsby 0;1i*gamab 0 epsbz ];

epsa=[epsax 0 -1i*gamaa;0 epsay 0;1i*gamaa 0 epsaz ];
else
epsb=[epsbx];

epsa=[epsax];
endif

line=getNewDataLine(fid);
numbs = str2num(line);
Na=numbs(1);


line=getNewDataLine(fid);
numbs = str2num(line);
transmit=numbs(1);
rotation=numbs(2);
numbs

line=getNewDataLine(fid);
numbs = str2num(line);

wn1=numbs(1);
wn2=numbs(2);
ndiv=numbs(3);

line=getNewDataLine(fid);
numbs = str2num(line);
nGx=numbs(1);
nGy=numbs(2);

line=getNewDataLine(fid);
numbs = str2num(line);
plotFT=numbs(1);
plotWave=numbs(2);

%=================

a=ax;

d1=0*a;
d2=ax/2-Rx+0*a;

t1=cputime;

dwn=(wn2-wn1)/ndiv;



cf=(plotWave==0);

for p=1:1*cf*ndiv+1
    
    p
    
    Fn(p)=wn1+dwn*(p-1);
    k1=2*pi*Fn(p)/a;
    
    
    [Ts Rs, Fs]=calculteFaraday(geometry,epsa,epsb,eps1,eps3,ax,ax,Rx,Ry,d1,d2,Na,nGx,nGy,k1,p,plotFT,plotWave,rotationGraphColor,theta,fi);
           
      if(real(Ts)>1) 
        Ts=1;
      end
      
    Tr(p)=real(Fs);
    Tt(p)=real(Ts);
    
    uu=Ts+Rs;

end

t2=cputime;

comptation_time=t2-t1;
comptation_time;
        
        Tr';
        Tt';
		
result=zeros(1*cf*ndiv+1,3);

  for p=1:1*cf*ndiv+1
  
  result(p,1)= Fn(p);
  result(p,2)= Tr(p);
  result(p,3)= Tt(p);
  end
  
  disp('Results:');
  disp('[wn *  Rotation * Transmitance ]');	
  disp(result);

  
if(rotation &&length(Tr)>1)
                figure(1)
             plot(Fn,Tr,rotationGraphColor);
 
             axis([wn1,wn2,-90,90]);
             hold on
end
            
if(transmit &&length(Tt)>1)
              figure(2)
               plot(Fn,Tt,transmitGraphColor);
                 axis([wn1,wn2,min(Tt),max(Tt)]);
                %axis([wn1,wn2,0,1]);
                 hold on
            
end

end


