function [Dc, Ac] =MeanValuefromCenter(A) %mean value from center of matrix. (A=2D matrix)
%%
[nx,ny]=size(A)
nxc=[nx/2];
nyc=ny/2;

for ii=1:nx
    for jj=1:ny
   D(ii,jj)=sqrt((nx-nxc-ii)^2+(ny-nyc-jj)^2);     
    end
end
   

Dl=D(:);
Al=A(:);
[Ds,I]=sort(Dl);
As=Al(I);
[Ac]=As(As~=0);
Dc=Ds(As~=0);


end
