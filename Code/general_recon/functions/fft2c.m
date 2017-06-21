function res = fft2c(x)

% res = fft2c(x)
% modified for more dimensions (eg coils)

imsize=size(x,1)*size(x,2);

for iter=1:size(x,3)
res(:,:,iter) = 1/sqrt(imsize)*fftshift(fft2(ifftshift(x(:,:,iter))));
end
% 