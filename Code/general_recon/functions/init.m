function param = init(P)

N=size(P.mask);
if P.xfmWeight>0
    %           XFM=Wavelet_SQKSP(size(mask),[1,2]);
    XFM=Wav_SPOT(N);   %15 times faster....
else
    XFM=IOP; %Identity
end

param.FT = []; % The measurement operator (undersmapled Fourier for example)
param.data = []; % measurements to reconstruct from

param.gradToll = 1e-30;	% step size tollerance stopping criterea (not used)

param.l1Smooth = 1e-15;	% smoothing parameter of L1 norm
param.pNorm = 1;  % type of norm to use (i.e. L1 L2 etc)

% line search parameters
param.lineSearchItnlim = 150;
param.lineSearchAlpha = 0.01;
% param.lineSearchAlpha=1e-5;
param.lineSearchBeta = 0.6;
param.lineSearchT0 = 1 ; % step size to start with

% initialize Parameters for reconstruction

param.XFM = XFM;
param.TV = TVOP;
param.TVWeight =P.TVWeight;     % TV penalty
param.TV2=TV2op;
param.TV2Weight=P.TGVfactor*param.TVWeight;

param.Itnlim = P.Itnlim;


if P.WeightedL2
    param.V=(P.MNSA.*P.mask);
    param.V=repmat((P.MNSA.*P.mask),[1 1 P.nc]);
else
    param.V=ones(size(P.MNSA,1),size(P.MNSA,2),P.nc);
end

if P.VNSAlambdaCorrection % to compare different experiments
    param.xfmWeight=P.xfmWeight*(mean(param.V(P.mask~=0)));
    param.TVWeight =param.TVWeight*(mean(param.V(P.mask~=0)));     % TV penalty
    param.TV2Weight=param.TV2Weight*(mean(param.V(P.mask~=0)));
else
    param.xfmWeight = P.xfmWeight;  % L1 wavelet penalty
end
param.Beta='PR_restart';
param.display=P.visualize_nlcg;
param.Debug=P.debug_nlcg;



