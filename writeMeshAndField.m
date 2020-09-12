function writeMeshAndField(nx,ny,nz, E,mode,Na)

factor=1;


Ne=nx*ny*nz;

Nn=Ne*8;
ie=0;

elVert=zeros(Ne,8);

for kz=1:nz
    for ky=1:ny
        for kx=1:nx
            ie=ie+1;
            for j=1:8
                elVert(ie,j)=(ie-1)*8+j;
            end
        end
    end
end




% 	Nn=(ny+1)*(nx+1)*(nz+1);
% 		Nx=nx+1;
% 		Nxy=(nx+1)*(ny+1);
% 		Ne=nx*ny*nz;
%
%
%     frx=1;
%     frxy=1;
%     ie=0;
%
%      elVert=zeros(Ne,8);

%     while(frxy*Nxy<Nn)
%    ie=ie+1;
%         elVert(ie,1)=frxy*Nxy+frx;
%           elVert(ie,2)=frxy*Nxy+1+frx;
%             elVert(ie,3)=frxy*Nxy+Nx+1+frx;
%               elVert(ie,4)=frxy*Nxy+Nx+frx;
%               elVert(ie,5)=(frxy-1)*Nxy+frx;
%                 elVert(ie,6)=(frxy-1)*Nxy+1+frx;
%                   elVert(ie,7)=(frxy-1)*Nxy+Nx+1+frx;
%                     elVert(ie,8)=(frxy-1)*Nxy+Nx+frx;
%
%
%
%
%         if (mod(frx+1,Nx)==0 &&mod(frx+Nx+1,Nxy)>0 )
%             frx=frx+2;
%         elseif (mod(frx+Nx+1,Nxy)==0 )
%
%             frxy=frxy+1;
%             frx=1;
%
%         else frx=frx+1;
%         end
%     end

if(mode>1)
    
    fid = fopen('C:\\Users\\Hassan\\\Documents\MATLAB\\myMesh.txt','wt');  % Note the 'wt' for writing in text mode
    
    
    
    fprintf(fid,'hexahedron\n');
    fprintf(fid,'//numberOfNodes\n');
    fprintf(fid,'%d\n',Nn);
    fprintf(fid,'//numberOfElements\n');
    fprintf(fid,'%d\n',Ne);
    fprintf(fid,'//numberOfRegs\n');
    fprintf(fid,'%d\n',1);
    fprintf(fid,'//factor\n');
    fprintf(fid,'%f\n',factor);
    
    for ie=1:Ne
        
        for k=1:7
            fprintf(fid,'%d,', elVert(ie,k));
        end
        fprintf(fid,'%d\n', elVert(ie,8));
    end
    
    
    
    
    y=linspace(0,Na,ny+1);
    
    %          x=linspace(0,.001,nx+1);
    %
    %          z=linspace(0,.001,nz+1);
    dx=.001;
    dy=y(2)-y(1);
    dz=.001;
    ddx=1.0/nx;
    D=(nx-1)*ddx/2;
    
    for kz=1:nz
        for ky=1:ny
            for kx=1:nx
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky),dz);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky),dz);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky)+dy,dz);
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky)+dy,dz);
                
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky),0);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky),0);
                fprintf(fid,'%f,%f,%f\n',dx+(kx-1)*ddx-D,y(ky)+dy,0);
                fprintf(fid,'%f,%f,%f\n',(kx-1)*ddx-D,y(ky)+dy,0);
            end
        end
    end
    
    fprintf(fid,'%d,%d,xxx',1,Ne);
    
    
    fclose(fid);
end



nL=size(E,2);
fid = fopen('myFields.txt','wt');  % Note the 'wt' for writing in text mode

fprintf(fid,'displacement\n');
fprintf(fid,'3\n');
fprintf(fid,'%d\n',Nn);

ie=0;

for ky=1:ny
    for kx=1:nx
        
        ie=ie+1;
        fprintf(fid,'%d\t%f\t%f\t%f\n',elVert(ie,1),E(kx,ky,1),E(kx,ky,2),E(kx,ky,3));
        
        fprintf(fid,'%d\t%f\t%f\t%f\n',elVert(ie,5),-E(kx,ky,1),-E(kx,ky,2),-E(kx,ky,3));
    end
end
fclose(fid);

end

