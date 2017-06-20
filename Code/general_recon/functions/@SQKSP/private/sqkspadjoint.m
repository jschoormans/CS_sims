function [kspace] = sqkspadjoint(ksp,dims,origsize);

if ndims(ksp)==3

S=size(ksp);
for ii=1:3; %only first 3 dimensions (spatial)
d=(S(ii)-origsize(ii))/2; 
c(ii,:)=[1+ceil(d),S(ii)-floor(d)];
end


kspace=ksp([c(1,1):c(1,2)],[c(2,1):c(2,2)],[c(3,1):c(3,2)]); 

else % 2D
    
S=size(ksp);
for ii=1:2; %only first 3 dimensions (spatial)
d=(S(ii)-origsize(ii))/2; 
c(ii,:)=[1+ceil(d),S(ii)-floor(d)];
end


kspace=ksp([c(1,1):c(1,2)],[c(2,1):c(2,2)]); 
    
    
end