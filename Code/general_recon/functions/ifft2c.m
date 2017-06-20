function res = ifft2c(x)
%rr
%
% res = ifft2c(x)
% 
% orthonormal centered 2D ifft
%
% (c) Michael Lustig 2005

% res = sqrt(length(x(:)))*fftshift(ifft2(ifftshift(x)));

res = sqrt(length(x(:)))*(ifft2(ifftshift(x))); % trying out stuff.
