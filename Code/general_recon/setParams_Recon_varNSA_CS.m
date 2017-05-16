    function P=setParams_Recon_varNSA_CS()
        P=struct;
        % [P.nx P.ny P.nz P.nc P.nNSA]=size(MR.Data); %could be placed elsewhere in
        % the code (if used)
        if ~isfield(P,'TVWeight');
            P.TVWeight = 1; end;	% Weight for TV penalty
        if ~isfield(P,'xfmWeight');
            P.xfmWeight = 0;end	% Weight for Transform L1 penalty
        if ~isfield(P,'Itnlim');
            P.Itnlim = 20;end;         % Number of iterations
%         if ~isfield(P,'reconslices');
%             P.reconslices = [1:P.nx];end; %could be placed elsewhere in
        % the code (if used)
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
        % if ~isfield(P,'resultsfolder')
        %     P.resultsfolder=uigetdir('select resultsfolder');end
        if ~ isfield(P,'savename')
            C=clock;
            P.savename=['recon-',date,'-',num2str(C(4)),'-',num2str(C(5))]; end
        if ~isfield(P,'sensemapsloop')
            P.sensemapsloop=0
        end
        
        if and((P.squareksp==false),(P.xfmWeight>0));
            error('wavelet enabled -- use squareksp!');end
    end