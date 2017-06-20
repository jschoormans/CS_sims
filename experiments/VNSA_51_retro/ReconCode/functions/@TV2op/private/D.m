function res = D(image)

%
% res = D(image)
%
% image = a 2D image
%
% This function computes the finite difference transform of the image
%
% Related functions:
%       adjD , invD 
%
%
% (c) Michael Lustig 2005

[sx,sy] = size(image);

Dx = image([2:end,end],:) - 2.*image + image([1,1:end-1],:);
Dy = image(:,[2:end,end]) - 2.*image + image(:,[1,1:end-1]);

%res = [sum(image(:))/sqrt(sx*sy); Dx(:);  Dy(:)]; 
res = cat(3,Dx,Dy);


