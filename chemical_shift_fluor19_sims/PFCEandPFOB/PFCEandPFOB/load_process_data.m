cd('L:\basic\divi\Projects\cosart\CS_simulations\chemical_shift_fluor19_sims\PFCEandPFOB\PFCEandPFOB')

%% LOAD BRUKER IMAGES
fid = fopen('\\amc.intra\users\J\jschoormans\home\Downloads\PFCEandPFOB\PFCEandPFOB\Readout1\2dseq');
R = fread(fid,'int16');
fclose(fid);

Im1=reshape(R,[128,128,128]);

fid = fopen('\\amc.intra\users\J\jschoormans\home\Downloads\PFCEandPFOB\PFCEandPFOB\Readout2\2dseq');
R = fread(fid,'int16');
fclose(fid);

Im2=reshape(R,[128,128,128]);

figure(1); 
subplot(121)
imshow(abs(squeeze(Im1(:,:,64))),[])
subplot(122)
imshow(abs(squeeze(Im2(:,:,64))),[])

%% load kspace readout 1

fid = fopen('\\amc.intra\users\J\jschoormans\home\Downloads\PFCEandPFOB\PFCEandPFOB\Readout1\fid');
R = fread(fid,'int32','l');
fclose(fid);

K1=reshape(R,[256,128,128]);
K1=permute(K1,[2 3 1]);
for i=1:128; for j=1:128; for k=1:2:256
KC1(i,j,ceil(k/2))=K1(i,j,k)+1i*K1(i,j,k+1);            
        end
    end
end

Recon1=ifftshift(ifftn(fftshift(KC1)));

imagine(abs(Recon1))

%% load kspace readout 2
fid = fopen('\\amc.intra\users\J\jschoormans\home\Downloads\PFCEandPFOB\PFCEandPFOB\Readout2\fid');
R = fread(fid,'int32','l');
fclose(fid);

K2=reshape(R,[256,128,128]);
K2=permute(K2,[2 3 1]);
for i=1:128; for j=1:128; for k=1:2:256
KC2(i,j,ceil(k/2))=K2(i,j,k)+1i*K2(i,j,k+1);            
        end
    end
end
KC2=permute(KC2,[3 2 1]);
Recon2=ifftshift(ifftn(fftshift(KC2)));

% imagine(abs(Recon2))

%% plot both kspaces next to each other to check if the orientations are correct!!!

figure(2); 
subplot(321)
imshow(abs(squeeze(Recon1(:,:,64))),[])

subplot(322)
imshow(abs(squeeze(Recon2(:,:,64))),[])

subplot(323)
imshow(abs(squeeze(Recon1(:,64,:))),[])

subplot(324)
imshow(abs(squeeze(Recon2(:,64,:))),[])

subplot(325)
imshow(abs(squeeze(Recon1(64,:,:))),[])

subplot(326)
imshow(abs(squeeze(Recon2(64,:,:))),[])

%% save kspaces
save('kspaces_readout1&2.mat','KC1','KC2')
