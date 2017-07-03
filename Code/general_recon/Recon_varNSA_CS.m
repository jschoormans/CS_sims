classdef Recon_varNSA_CS < MRecon
    % J SCHOORMANS MAY 2017 
    % IMPLEMENTS ALL STEPS OF VARIABLE NSA CS RECONSTRUCTION
    
    
properties
% No additional properties needed
P=setParams_Recon_varNSA_CS(); % change this

end
methods
    function MR = Recon_varNSA_CS( filename )
        MR = MR@MRecon(filename);
        MR.P.resize=true;
        MR.P.resize_size=3;
        MR.P.fixsliceintensities=false; 

        MR.P.visualize_nlcg=0;
        MR.P.debug_nlcg=0;

        MR.P.parallel=false;
        MR.P.VNSAlambdaCorrection=1;
        MR.P.WeightedL2=1;
        MR.P.TVWeight=1e-5;
        MR.P.xfmWeight=2e-5;
        MR.P.Itnlim=25;
        MR.P.outeriter=1;
        MR.P.TGVfactor=0;
        MR.Parameter.Recon.RemovePOversampling='No';
        MR.P.prewhiten=true; 
        MR.Parameter.Recon.ImmediateAveraging='No' ;
        MR.Parameter.Recon.ArrayCompression='No'; %MUST BE NO, for prewhitening 
        assert(~(strcmp(MR.Parameter.Recon.ArrayCompression,'Yes')&&MR.P.cc),'Impossible, MR.Parameter.Recon.ArrayCompression should be : NO');
        
        MR.P.cc_number=8;
        MR.P.cc=true;
        MR.P.cctype='E'; %E,S,A,G (to test!)
        MR.P.Scaling=true;
    
        MR.P.VNorm=1; %power of weighting matrix (number of NSA)^p; 
    end
    % Overload (overwrite) the existing Perform function of MRecon
        function Perform( MR )
           MR.Perform1;
           MR.ReconCS;
           MR.fixsliceintensities;
           MR.resizerecon;

        end

    
    function Perform1( MR )   
        % Produce k-space Data (using existing MRecon functions)
        MRn=MR.Copy; MRn.Parameter.Parameter2Read.typ=5;
        MRn.ReadData;
        MR.P.eta=MRn.Data;

        MR.Parameter.Parameter2Read.typ=1;
        MR.ReadData;
        MR.DcOffsetCorrection;
        MR.PDACorrection;
        MR.RandomPhaseCorrection;
        MR.MeasPhaseCorrection;
        MR.addAveragestoLabels
        MR.SortData;
        MR.Data=MR.Data;
        [MR.P.mask,MR.P.MNSA,MR.P.pdf]=makemask(MR);
        MR.Average;      
        MR.preWhiten;

        checkerboard=create_checkerboard([1,size(MR.Data,2),size(MR.Data,3)]);
        MR.Data=bsxfun(@times,checkerboard,MR.Data);    %undo checkerboard    
        
        if MR.P.cc        % coil compression 
        MR.Data=bart(['cc -',MR.P.cctype,' -p ',num2str(MR.P.cc_number)],MR.Data);
        end
        
        MR.P.sensemaps=estsensemaps(MR);        % calculate sense maps
        MR.Data=fftshift(ifft(ifftshift(MR.Data,1),[],1),1);
        MR.ScaleKSP;
        
        %remove oversampling from sense maps and MR.Data;
        nsl=size(MR.Data,1)/MR.Parameter.Encoding.KxOversampling;
        cropslices=floor((size(MR.Data,1)/2)-nsl/2):floor((size(MR.Data,1)/2)-nsl/2)+nsl-1;
        MR.P.sensemaps=MR.P.sensemaps(cropslices,:,:,:);
        MR.Data=MR.Data(cropslices,:,:,:);
        
%         MR.Data=MR.Data./max(abs(MR.Data(:))); %normalize
        
        % add other parameters in P
        [MR.P.nx MR.P.ny MR.P.nz MR.P.nc MR.P.nNSA]=size(MR.Data); %could be placed elsewhere in
        if ~isfield(MR.P,'reconslices'); MR.P.reconslices=[1:MR.P.nx];        end
        disp('Perform 1 finished!')
        
    end
    function ReconCS(MR)
        recon=zeros(length(MR.P.reconslices),size(MR.Data,2),size(MR.Data,3)); %pre-allocate

        if MR.P.parallel==true
            
            tic
            parfor sl_iter=1:length(MR.P.reconslices) %loop over slices
                sl=MR.P.reconslices(sl_iter);
                recondata=MR.Data(sl,:,:,:); % data of one slice to be used in recon
                sensemapsslice=squeeze(MR.P.sensemaps(sl,:,:,:));
                param=init(MR.P);
                param_sl=setReconParams(MR,recondata,sensemapsslice,param);
                recon(sl_iter,:,:)=runCS(MR,param_sl,MR.P);
            end
            disp('CS finished!')
            toc
        else
            param=init(MR.P);
            for sl_iter=1:length(MR.P.reconslices)  %loop over slices
                sl=MR.P.reconslices(sl_iter);
                tic;
                recondata=MR.Data(sl,:,:,:); % data of one slice to be used in recon
                sensemapsslice=squeeze(MR.P.sensemaps(sl,:,:,:));
                param_sl=setReconParams(MR,recondata,sensemapsslice,param);
                recon(sl_iter,:,:)=runCS(MR,param_sl,MR.P);
                time_iter=toc;
                fprintf('slice: %u of %u | delta t: %f2 \n',sl,length(MR.P.reconslices),time_iter)
            end
        end
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
        for iline=numbers(MR.Parameter.Labels.Index.typ==1);
            chan=MR.Parameter.Labels.Index.chan(iline);
            xc=all_klines(iline,1)+1+abs(mk(1));
            yx=all_klines(iline,2)+1+abs(mk(2));
            nsa=(V(xc,yx,chan));
            V(xc,yx,chan)=V(xc,yx,chan)+1; % is actually MNSA --> easier to keep other function for now
            nsa_number(iline) = nsa;
        end
        MR.Parameter.Labels.Index.aver=(nsa_number).';
        MR.Parameter.Labels.Index.dyn=zeros(size(MR.Parameter.Labels.Index.dyn)); %remove dynamics
    end % GIVES BUGS NOW???
    
    function [mask,MNSA,pdf]=makemask(MR)
        K=MR.Data;
        % 1) make mask
        disp('make mask and setting up data...')
        Ks=squeeze(K(round(size(K,1)/2),:,:,1,:)); %k-space for one channel and one slice
        fullmask=Ks~=0;             %find mask used for scan (nx*ny*NSA)
        MNSA=sum(fullmask,3);       %find NSA for all k-points in mask
        
        % 4) calculate 2D mask (??) and PDF
        mask=double(MNSA>0);        %2d mask (no NSA dimension)
        pdf=estPDF(mask);       %pdf is used for first guess; should be fixed!
    end
    
    function pdf=estPDF(mask)
        h=1/49*ones(7);
        pdf=conv2(mask,h,'same');
        pdf=pdf+eps;
    end
    
    function sensemaps=estsensemaps(MR)
        disp('Estimating sense maps...'); tic;
        sensemaps=bart('ecalib -m1 -S -I -c0.01',MR.Data);toc;
        sensemaps=conj(sensemaps);
        
        if MR.P.Scaling        % scale sense maps
            disp('Scaling sense maps...')
            [nx,ny,nz,nc]=size(sensemaps); 
            sensmag=reshape(vecnorm(reshape(sensemaps,[],nc).'),[nx ny nz]);
            sensemaps=bsxfun(@rdivide,sensemaps,sensmag);
        end
        
    end
    
    function param=setReconParams(MR,recondata,sensemaps,param)
%         param = init(MR.P)
        
        param.data=squeeze(recondata);
        N=size(MR.P.mask);
        % low res phase estimate        
        FT1 = MCp2DFT(MR.P.mask, N,sensemaps, 1, 2);       
        I=FT1'*param.data; 
        ph=angle(I);
        ph=exp(1i*ph);
        FT = MCp2DFT(MR.P.mask, N,sensemaps, ph, 2);
        param.FT = FT;

    end
    
    function recon=runCS(MR,param,P)
        res=param.XFM*(param.FT'*(param.data./repmat(P.pdf,[1 1 P.nc])));
%           res=param.XFM*(param.FT'*(param.data));

        for n=1:P.outeriter
            res = MCfnlCg_test(res,param);
        end
        recon = param.XFM'*res;
%         recon=recon./max(recon(:));
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
    
    function preWhiten(MR)
        if MR.P.prewhiten
            % after sortdata --> (?) before averaging (?)
            disp('pre-whitening data...')
            sizeMRDATA=size(MR.Data);
            Ncoils=size(MR.Data,4);
            Nsamples=numel(MR.Data)/Ncoils;
            data=reshape(MR.Data,[Nsamples,Ncoils]);
            data=data.';
            
            psi = (1/(Nsamples-1))*(MR.P.eta' * MR.P.eta);
            L = chol(psi,'lower');
            L_inv = inv(L);
            data_corr = L_inv * data;
            data_corr=data_corr.';
            MR.Data=reshape(data_corr,sizeMRDATA);
            disp('done')
        else
        end

    end
    
    function ScaleKSP(MR)
        fprintf('\n Scaling KSP...')
        if MR.P.Scaling
            tic % maybe takes too long with the loop.
            ksptemp1=permute(MR.Data,[2 3 4 1]);
            ksptemp=reshape(MR.Data,size(ksptemp1,1),size(ksptemp1,2),[]);
            ksptemp=ifft2c(ksptemp); %dont really care about ffshift for this
            ksptemp=reshape(ksptemp,size(ksptemp1));
            ksptemp=permute(ksptemp,[4 1 2 3]);
            toc 
            tmp = dimnorm(ksptemp, 4);
            tmpnorm2 = sort(tmp(:), 'ascend');
            p100 = tmpnorm2(end);
            p90 = tmpnorm2(round(.9 * length(tmpnorm2)));
            p50 = tmpnorm2(round(.5 * length(tmpnorm2)));
            if (p100 - p90) < 2 * (p90 - p50)
                scaling = p90;
            else
                scaling = p100;
            end
            MR.P.scaling=scaling; % save for later use
            fprintf('\nScaling: %f\n\n', scaling);
            
            MR.Data = MR.Data ./ scaling;
        end
    end
    
end
end