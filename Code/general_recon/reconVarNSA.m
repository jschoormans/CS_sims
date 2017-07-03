function P=reconVarNSA(K,P)
%   JASPER SCHOORMANS 25-10-2016
%   RECONSTRUCTION OF VARIABLE NSA CS MEASUREMENTS
%   INPUT: K a k-space matrix [nx ny nz nc nNSA]

P=setParams(K,P);
% K=FFTmeas(K,P);
[data,P.mask,P.MNSA,P.pdf]=makemask(K,P);
for sl=P.reconslices %loop over slices
recondata=data(sl,:,:,:); % data of one slice to be used in recon 
param=setReconParams(recondata,P.MNSA,P.mask,P.pdf,P.sensemaps(:,:,:),P);
recon(:,:,sl)=runCS(param,P);
end

P.recon=recon;
% saverecon(recon,P)



end

function P=setParams(K,P)
[P.nx P.ny P.nz P.nc P.nNSA]=size(K);
if ~isfield(P,'TVWeight');
    P.TVWeight = 0.000; end;	% Weight for TV penalty
if ~isfield(P,'xfmWeight');
    P.xfmWeight = 0.1;end	% Weight for Transform L1 penalty
if ~isfield(P,'Itnlim');
    P.Itnlim = 20;end;         % Number of iterations
if ~isfield(P,'reconslices');
    P.reconslices = [1:P.nx];end; 
if ~isfield(P,'squareksp')%makes ksp bigger until first 2^n
    P.squareksp=true;end
if ~isfield(P,'outeriter')%
    P.outeriter=6;end
if ~isfield(P,'break')%
    P.break=false;end
if ~isfield(P,'noNSAcorr')%
    P.noNSAcorr=false;end
if ~isfield(P,'TGVfactor')%
    P.TGVfactor=2;end
if ~isfield(P,'parfor')%
    P.parfor=0;end
if ~isfield(P,'resultsfolder')
    P.resultsfolder=uigetdir('select resultsfolder');end
if ~ isfield(P,'savename')
    C=clock;
    P.savename=['recon-',date,'-',num2str(C(4)),'-',num2str(C(5))]; 
end
if ~isfield(P,'sensemapsloop')
    P.sensemapsloop=0
end


end


function K=FFTmeas(K,P)
% FFT IN MEASUREMENT DIRECTIONS+normalization
tic; disp('FFT in measurement direction')
K=ifftshift(ifft(K,[],1),1);;
K=K./max(K(:)); %normalize kspace
%K=squeeze(K); %WHY THIS???
disp('if errors with real data uncomment squeeze' )
toc;

end

function [data,mask,MNSA,pdf]=makemask(K,P)

% 1) make mask
tic; disp('make mask and setting up data...')
Ks=squeeze(K(round(size(K,1)/2),:,:,1,:)); %k-space for one channel and one slice
fullmask=Ks~=0;             %find mask used for scan (nx*ny*NSA)
MNSA=sum(fullmask,3);       %find NSA for all k-points in mask 

% 2) average over NSA 
data=sum(K,5);              %sum of data over NSA
data=data./permute(...
    repmat(MNSA,[1 1 P.nx P.nc]),[3 1 2 4]);            %mean of data over NSA 2
data(isnan(data))=0;        %clear up NaN values (due to /0)

% 3) make square KSPACE
if P.squareksp==true
data=squareksp(data,[2 3]);       %make k-spsace square and size a 2^n (bit buggy, not needed for nx now)
MNSA=squareksp(MNSA);       %make MNSA square and size of 2^n
end

% 4) calculate 2D mask (??) and PDF
mask=double(MNSA>0);        %2d mask (no NSA dimension)
pdf=estPDF(mask);       %pdf is used for first guess; should be fixed!

disp(['acceleration: ',num2str(sum(mask(:))/(P.ny*P.nz))])
toc
end

function pdf=estPDF(mask)
h=1/49*ones(7);
pdf=conv2(mask,h,'same');
pdf=pdf+eps;
end

function sensemaps=estsensemaps(recondata,P)
if ~isfield(P,'sensemapsprovided') | P.sensemapsprovided==0;
    disp('Estimating sense maps');tic
    
    if P.sensemapsloop==1
    if true
        recon=bart('fftmod 7',bart('fft 7',bart('fftmod 7',recondata)));
        ksp=bart('fftmod -i 7',bart('fft -i 7',bart('fftmod -i 7',recon)));
        sensemaps=bart('ecalib -r15 -m1',ksp);
        sensemaps=fftshift(fftshift(sensemaps,2),3);
        sensemaps=(sensemaps);
        %FFTMOD???
    else
        %     sensemaps=bart('ecalib -r15 -m1',bart('fftmod -i 7',recondata));
        sensemaps=bart('ecalib -r15 -m1',recondata);

    end
    else
        recon=bart('fftmod 7',bart('fft 7',bart('fftmod 7',recondata)));
        ksp=bart('fftmod -i 7',bart('fft -i 7',bart('fftmod -i 7',recon)));
        sensemaps=bart('ecalib -r15 -m1',ksp);
    end
    toc
else
    disp('using provided sense maps')
    sensemaps=P.sensemaps;
end
if ndims(squeeze(sensemaps))==3 % for 2D
    %     sensemaps=permute(sensemaps,[4 1 2 3])
    sensemaps=fftshift(sensemaps,3);
    sensemaps=fftshift(sensemaps,2);
    
end

if P.squareksp==true; 
sensemaps=squareksp(sensemaps);     
end

end

function param=setReconParams(recondata,MNSA,mask,pdf,sensemaps,P)
tic; disp('setting l1-recon parameters');
N=size(mask); 
param = init(P);

sensemaps=squeeze(sensemaps); 
param.data=squeeze(recondata);

N=size(P.mask);
if P.xfmWeight>0
    %           XFM=Wavelet_SQKSP(size(mask),[1,2]);
    XFM=Wav_SPOT(N);   %15 times faster....
else
    XFM=IOP; %Identity
end


FT1 = MCp2DFT(mask, N,squeeze(conj(sensemaps)), 1, 2);
angle_estimate=FT1'*param.data;

ph = angle(angle_estimate);
ph=exp(i*ph);
FT = MCp2DFT(mask, N,squeeze(conj(sensemaps)), ph, 2);
% imshow(angle(FT'*(param.data)),[])
% initialize Parameters for reconstruction
param.XFM = XFM; %easiest removal is to replace with empty operator???

param.TV = TVOP;
param.TVWeight =P.TVWeight;     % TV penalty
param.TV2=TV2op;
param.TV2Weight=P.TGVfactor*param.TVWeight;

param.Itnlim = P.Itnlim;
param.lineSearchItnlim=100;
param.Debug=0;
param.lineSearchAlpha=1e-5;

if P.noNSAcorr
    param.V=ones(size(MNSA,1),size(MNSA,2),P.nc);
    param.xfmWeight = P.xfmWeight;  % L1 wavelet penalty
    
else
    param.V=(MNSA.*mask);
    param.xfmWeight=P.xfmWeight*(mean(param.V(mask~=0)));
    param.TVWeight =param.TVWeight*(mean(param.V(mask~=0)));     % TV penalty
    param.TV2Weight=param.TV2Weight*(mean(param.V(mask~=0)))
    param.V=repmat((MNSA.*mask),[1 1 P.nc]);
end
param.FT = FT;
param.Beta='PR_restart';
param.display=1;
param.display=P.visualize_nlcg;
param.Debug=P.debug_nlcg;
toc;
end

function recon=runCS(param,P)
res=param.XFM*(param.FT'*(param.data./repmat(P.pdf,[1 1 P.nc])));
for n=1:P.outeriter
    res = MCfnlCg_test(res,param);
end
recon = param.XFM'*res;
recon=recon./max(recon(:));
end

function saverecon(recon,P)
save recons (TODO!!)
disp('saving recons...')
cd(P.resultsfolder)

ni=make_nii(abs(recon))
save_nii(ni,[P.savename,'.nii'])

% save(P.savename,'recon')
% disp('saving images...')
% A=figure(5); imshow(abs(recon));axis off;
% export_fig(P.savename,'-native')
% B=figure(6); imshow(abs(P.MNSA));axis off;colormap jet; colorbar
% export_fig([P.savename,'mask'],'-native')

diary([P.savename,'settings.txt']);
disp('saving settings...')
P
disp('Finished!'); diary off
end

function ph = phase_estimate(recondata,sensemaps,P)
% low res phase estimate 
I=bart(['resize -c 1 ',num2str(floor(size(recondata,2)/4)),' 2 ',num2str(floor(size(recondata,3)/4))],recondata);
I=bart(['resize -c 1 ',num2str(size(recondata,2)),' 2 ',num2str(size(recondata,2))],I);
I=bart('fft -i 7',I);
I=ifftshift(ifftshift(I,2),3);
I=bart('fmac -C',I,sensemaps);
I=sum(I,4);
I=squeeze(I);
ph=angle(I);
end