classdef Recon_varNSA_CS < MRecon
    % J SCHOORMANS MAY 2017 
    % IMPLEMENTS ALL STEPS OF VARIABLE NSA CS RECONSTRUCTION
    
    
properties
% No additional properties needed
P=setParams_Recon_varNSA_CS();

end
methods
    function MR = Recon_varNSA_CS( filename )
        % Create an MRecon object MR upon creation of My_Recon
        MR = MR@MRecon(filename);
    end
    % Overload (overwrite) the existing Perform function of MRecon
        function Perform( MR )
           MR.Perform1
           MR.ReconCS
           MR.fixsliceintensities
           MR.resizerecon;

        end

    
    function Perform1( MR )
        MR.Parameter.Parameter2Read.typ = 1;
        MR.Parameter.Recon.ImmediateAveraging='No' 
               
        % Produce k-space Data (using existing MRecon functions)
        MR.ReadData;
        MR.DcOffsetCorrection;
        MR.PDACorrection;
        MR.RandomPhaseCorrection;
        MR.MeasPhaseCorrection;
        MR.addAveragestoLabels
        
        disp('sortdata')
        MR.SortData;
        
        % add other parameters in P
        [MR.P.nx MR.P.ny MR.P.nz MR.P.nc MR.P.nNSA]=size(MR.Data); %could be placed elsewhere in
        if ~isfield(MR.P,'resize'); MR.P.resize=true; end; 
        if ~isfield(MR.P,'resize_size'); MR.P.resize_size=3; end; 
        if ~isfield(MR.P,'MR.P.fixsliceintensities'); MR.P.MR.P.fixsliceintensities=true; end
        % make MNSA mask (could we do this before SortData)?
        disp('make mask...')
        [MR.P.mask,MR.P.MNSA,MR.P.pdf]=makemask(MR);
        
        % do averaging
        MR.Average;
        
        % make square kspace (if this is one of the options?)
        
        
        % calculate sense maps
        sensemaps=estsensemaps(MR);
        MR.P.sensemaps=sensemaps;
        MR.K2IM
        
        MR.Data=MR.Data./max(abs(MR.Data(:))); %normalize
        disp('Perform 1 finished!')
%         if MR.P.squareksp==true
%             MR.Data=squareksp(MR.Data,[2 3]);       %make k-spsace square and size a 2^n (bit buggy, not needed for nx now)
%             MR.P.sensemaps=squareksp(MR.P.sensemaps,[2,3]);
%             MR.P.mask=squareksp(MR.P.mask);
%             MR.P.MNSA=squareksp(MR.P.MNSA);
%             MR.P.pdf=squareksp(MR.P.pdf);
%         end
        
        
    end
    function ReconCS(MR)
        
        if isempty('MR.P.reconslices')
        MR.P.reconslices=[1:MR.P.nx];
        end
        
        % for all slices: set reconparams
        for sl=MR.P.reconslices %loop over slices
            recondata=MR.Data(sl,:,:,:); % data of one slice to be used in recon
            sensemapsslice=squeeze(MR.P.sensemaps(sl,:,:,:));
            param=setReconParams(MR,recondata,MR.P.MNSA,MR.P.mask,MR.P.pdf,sensemapsslice,MR.P);
            recon(:,:,sl)=runCS(MR,param,MR.P);
        end
%         MR.Data=recon;  % UNCOMMENT WHEN CODE IS TRULY FINISHED 
        MR.P.Recon=recon;
    end
    
    function addAveragestoLabels(MR)
        
        all_klines=double([MR.Parameter.Labels.Index.ky MR.Parameter.Labels.Index.kz]);
        all_klines(MR.Parameter.Labels.Index.typ~=1)=NaN;
        numbers=1:length(MR.Parameter.Labels.Index.typ);
        
        mm=length(numbers);
        mk=min(all_klines);
        mak=max(all_klines);
        
        nchans=max(MR.Parameter.Labels.Index.chan);
        
        V=zeros(mak(1)+abs(mk(1))+1,mak(2)+abs(mk(2))+1,nchans);
        tic
        for iline=numbers(MR.Parameter.Labels.Index.typ==1);
            
            chan=MR.Parameter.Labels.Index.chan(iline);
            
            xc=all_klines(iline,1)+1+abs(mk(1));
            yx=all_klines(iline,2)+1+abs(mk(2));
            
            nsa=(V(xc,yx,chan));
            V(xc,yx,chan)=V(xc,yx,chan)+1; % is actually MNSA --> easier to keep other function for now
            
            nsa_number(iline) = nsa;
        end
        toc
        MR.Parameter.Labels.Index.aver=(nsa_number).';
        MR.Parameter.Labels.Index.dyn=zeros(size(MR.Parameter.Labels.Index.dyn)); %remove dynamics
    end % GIVES BUGS NOW???
    
    function [mask,MNSA,pdf]=makemask(MR)
        K=MR.Data;
        % 1) make mask
        tic; disp('make mask and setting up data...')
        Ks=squeeze(K(round(size(K,1)/2),:,:,1,:)); %k-space for one channel and one slice
        fullmask=Ks~=0;             %find mask used for scan (nx*ny*NSA)
        MNSA=sum(fullmask,3);       %find NSA for all k-points in mask
        
        % 4) calculate 2D mask (??) and PDF
        mask=double(MNSA>0);        %2d mask (no NSA dimension)
        pdf=estPDF(mask);       %pdf is used for first guess; should be fixed!
        toc
    end
    
    function pdf=estPDF(mask)
        h=1/49*ones(7);
        pdf=conv2(mask,h,'same');
        pdf=pdf+eps;
    end
    
    function sensemaps=estsensemaps(MR)
        disp('estimating sense maps...'); tic;
        sensemaps=bart('ecalib -m1 -S -I -c0.01',MR.Data);toc;
        sensemaps=fftshift(fftshift(sensemaps,2),3);
        sensemaps=conj(sensemaps);
    end
    
    function param=setReconParams(MR,recondata,MNSA,mask,pdf,sensemaps,P)
        tic; disp('setting l1-recon parameters');
        addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/sparseMRI_v0.2'));
        addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/Wavelab850'));
        addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/spot-master'))

        
        N=size(mask);
        
        %generate transform operator
        
        if P.xfmWeight>0
            XFM=Wavelet_SQKSP(size(mask),[1,2]);
        else
            XFM=IOP; %Identity 
        end
        param = init;
        param.data=squeeze(recondata);

        % low res phase estimate        
        FT1 = MCp2DFT(mask, N,sensemaps, 1, 2);       

        I=FT1'*param.data; 
        ph=angle(I);
        ph=exp(1i*ph);
        FT = MCp2DFT(mask, N,sensemaps, ph, 2); param.FT = FT;
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
        param.Beta='PR_restart';
        param.display=1;
        
        toc;
    end
    
    function recon=runCS(MR,param,P)
        res=param.XFM*(param.FT'*(param.data./repmat(P.pdf,[1 1 P.nc])));
        for n=1:P.outeriter
            res = MCfnlCg_test(res,param);
        end
        recon = param.XFM'*res;
        recon=recon./max(recon(:));
    end
    
    
    function ph = phase_estimate(recondata,sensemaps)
        % low res phase estimate
        I=bart('resize -c 1 100 2 100',recondata);  
        I=bart('resize -c 1 1024 2 1024',I);
        
        I=bart('fft -i 7',I);
        I=ifftshift(ifftshift(I,2),3);
        I=bart('fmac -C',I,sensemaps);
        I=sum(I,4);
        I=squeeze(I);
        ph=angle(I);
    end
    
    function fixsliceintensities(MR)
        if MR.P.fixsliceintensities==true
        disp('smoothing slice intensities')
        recon=MR.P.Recon;
        [nx ny nz]=size(recon);
        for i=1:nz
            scaling(i)=mean(mean(squeeze(abs(recon(:,:,i)))));
        end
        scalings=smooth(scaling,15);
        S=scalings./(scaling+eps).';
        for i=1:nz
            recon(:,:,i)=recon(:,:,i).*S(i);
        end
        MR.P.Recon=recon;
        end
    end
    
    
    function resizerecon(MR)
        % interpolates result for better image quality by zerofilling in
        % kspace
        disp('interpolating...')
        recon=MR.P.Recon;
        if MR.P.resize==true
        
        x1=round(size(recon,1)*MR.P.resize_size);
        y1=round(size(recon,2)*MR.P.resize_size);
        z1=round(size(recon,3)*MR.P.resize_size);
        krr=bart('fft 7',recon);
        krr=bart(['resize 0 ',num2str(x1),' 1 ',num2str(y1),' 2 ',num2str(z1)], krr);
        krr=bart('fft -i 7',krr);
        MR.P.Recon=krr;      
        end
    end
end
end