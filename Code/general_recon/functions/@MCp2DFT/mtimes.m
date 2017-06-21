function res = mtimes(a,b)
% res = mtimes(FT, x)
%


if a.adjoint %FROM MULTICOIL KSPACE TO COMBINED IMAGE SPACE
    for ch=1:1:a.nchans %will size work tho?
    bb = reshape(b(:,:,ch),a.dataSize(1),a.dataSize(2));
%     resi = zpad(bb.*a.mask,a.imSize(1),a.imSize(2));
    resi=bb.*a.mask; 
    resi = ifft2c(resi);
    res(:,:,ch) = resi.*conj(a.ph);
    end
    
    res=sum(res.*a.sensmaps,3)./sum(abs(a.sensmaps+eps).^2,3);
%     res=res.*conj(a.ph)
else
    b = reshape(b,a.dataSize(1),a.dataSize(2));
    b=b.*(a.ph); %18-5-2017 added to maintain good phase for many iterations

    resi=bsxfun(@times,b,conj(a.sensmaps));
    res=bsxfun(@times,fft2c(resi),a.mask);
end

if size(b,2) == 1
    res = res(:);
end


