

MR=Recon_varNSA_CS(strcat(folder,files(iter).name));
MR.Parameter.Recon.ArrayCompression='No'; %coil compression
MR.Parameter.Recon.ACNrVirtualChannels=4;
MR.ReadData;
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
%%
MR.SortData;


%% calculate noise covariance of data
data=squeeze(MR.Data{1});
size(data)
data1=squeeze(data(:,:,:,:));
data1=reshape(data1,[576*256*85,12]);
size(data1);
c1=cov(data1);
size(c1)

%% calculate noise covariance of noise data
eta=squeeze(MR.Data{5});
size(eta)
c2=cov(eta);
size(c2);

figure(1); 
subplot(411);imshow(c1,[]); colormap('jet'); 
title('covariance matrix of data');
subplot(412);imshow(c2,[]); colormap('jet');
title('covariance matrix of noise')


%% 
Nsamples=576*256*85;
Ncoils=12;
psi = (1/(Nsamples-1))*(eta' * eta);

subplot(413);imshow(psi,[]); colormap('jet');
title('Psi')
%%


L = chol(psi,'lower');
L_inv = inv(L);

subplot(414);imshow(L_inv,[]); colormap('jet');
title('L_inv')
%%
data_corr = L_inv * data1.';
