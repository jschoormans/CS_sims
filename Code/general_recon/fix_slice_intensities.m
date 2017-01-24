% intensity correction for nifftii

cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/results_VNSA2')
D=dir('*.nii')

for ii=1:7
    close all
    cd('/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/results_VNSA2')

N=load_nii(D(ii).name)
%%
[nx ny nz]=size(N.img);

for i=1:nz
    scaling(i)=mean(mean(N.img(:,:,i)));
end
scalings=smooth(scaling,15)

S=scalings./scaling.';
figure(1);
hold on 
plot(scaling); plot(scalings)
hold off

clear I
for i=1:nz
    i
I(:,:,i)=N.img(:,:,i).*S(i);
end

figure(2);
subplot(211)
imshow(abs(squeeze(N.img(150,:,:))),[])

subplot(212)
imshow(abs(squeeze(I(150,:,:))),[])

filename=[D(ii).name(1:26),'_fixed.nii']
nii=make_nii(I)
save_nii(nii,filename)
pause(5)
end
