function Mfull = genEllipse(ny,nz)
%generates an elliptical k-space sampling pattern

Mfull=zeros(ny,nz);
for i=1:ny;
    for j=1:nz;
Mfull(i,j)=(i-0.5-(ny/2))^2+((ny/nz)*(j-0.5-(nz/2)))^2<(ny/2)^2 ;
end
end
