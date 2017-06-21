function res = ifft2c(x)
% res = ifft2c(x)
% modified for more dimensions (eg coils)
imsize=size(x,1)*size(x,2);

for iter=1:size(x,3)
res(:,:,iter) = sqrt(size(x,1)*size(x,2))*fftshift(ifft2(ifftshift(x(:,:,iter))));
end
