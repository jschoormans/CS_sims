[PSF,OTF]=PSFinF(128,PFOB,PFOB_alpha,1);
OTFl=repmat(OTF,[128 1]);
%
D=opDiag(OTFl(:))
deconv_image1=opInverse(F)*D*(data(1:128^2));
% deconv_image1=opInverse(F)*(data(1:128^2));

figure(52); imshow(abs(rr1(deconv_image1)),[])

%%

[PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(282.56877,44642/128) %126.4??
% [PFCE,PFCE_alpha,PFOB,PFOB_alpha]=calcspectra_BW(282.56877,44442/128) %126.4??

PFOB_alpha=PFOB_alpha%./(sum(PFOB_alpha))
[PSF,OTF]=PSFinF(128,-PFOB,PFOB_alpha,1);
OTFl=repmat(OTF,[128 1]);
OTFr=OTFl.';
D=opDiag(OTFl(:));
D2=opDiag(OTFr(:));

Spectrum2=[PFCE_alpha]; 
B=opConvolve(128,128,Spectrum2,[0 0],'cyclic') 
B2=opConvolve(128,128,Spectrum2.',[0 0],'cyclic')
FS=opConvolve(128,128,[1],[0 0],'cyclic')

M1=[D*F*FS,F*B*FS];
M2=[D2*F*FS,F*B2*FS];
M=[M1;M2] %measurement operator 

figure(99); 
subplot(211);imshow(abs(rr1(opInverse(D*F)*data(1:128^2,:))),[0 1e-2])
subplot(212);imshow(abs(rr(opInverse(M)*data)),[0 1e-2])

%%
CGK=zeros(size(data));
for outeriter=1:5
    CGK=nl_conjgrad_fluor(M,data,(CGK),20,zeros(size(data)),1e-4,128,256);
end
figure(2); subplot(224); imshow(rr(abs(CGK)),[]); axis off; title('CG-k space')