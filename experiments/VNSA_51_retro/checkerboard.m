close all

FT = MCp2DFT(mask, N,squeeze(ones(size((sensemaps)))), ph, 2);

if false %2d checkerboard
checkerboard=double(invhilb(1024) > 0);
checkerboard=checkerboard-double(checkerboard==0);
data= repmat(checkerboard,[1 1 8]).* param.data;
elseif false% 1d checkerboard
    rr=mod([1:1024],2); rr=double(rr)-double(rr==0);
    checkerboard=repmat(rr,[1024 1]);
    data= repmat(checkerboard,[1 1 8]).* param.data;
else
data=param.data;    
end

% data=param.data;

figure(1)
imshow(abs(FT'*data(:,:,:)),[])
figure(2)
imshow(abs(data(:,:,1)),[])
figure(3)
imshow(angle(data(:,:,1)),[-pi, pi])

figure(4)
imshow(abs(sensemaps(:,:,1)),[])
