function x=MultiCoil_CG(x0,params)

%   MODIFICATION OF CONJUGATE GRADIENT ALGO BY M LUSTIG
%   USING BART TOOLBOX
%   SIMPLIFIED VERSION
%   non-linear conjugate gradient algorithm that uses variable density
%   sampling
%   dependencies: BART toolbox

x = x0;

% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;
params.mask=params.data~=0;

g0 = wGradient(x,params); %first gradient at x0 
dx = -g0;                 %dx is minus g ??


% iterations
while(1)
    lsiter = 0;

    % pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, 0, params); %objective at original point (line length 0)
	t = t0;
    [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params); %objective along line 

    
    
    
    while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter)  %stopping critertion
        lsiter = lsiter + 1;
        t = t * beta;
        [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
    end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2       %line search took long; new t0 is made shorter
		t0 = t0 * beta
	end 
	
	if lsiter<1         %line search converged too fast: new t0 is longer
		t0 = t0 / beta
	end

	x = (x + t*dx);     %new image x is old x plus t*dx

	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d, ErrObj: %f', k,f1,RMSerr,lsiter,ERRobj));
    figure(100); imshow(abs(squeeze(x(1,:,:))),[])
	%---------------------------------------------------------------
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    bk=0 %steepest descent search
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
	end

end

return;
end

function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap

Im=bart('cdf97 -i 6',x);
Imdx=bart('cdf97 -i 6',dx);

FTXFMtx = bart('fakeksp',bart('fftmod -i 6',Im),ones(size(params.data)),params.sensemaps);        %multi-coil kspace of X
FTXFMtx = bart('fftmod -i 6',FTXFMtx);
FTXFMtx=FTXFMtx.*params.mask;
FTXFMtdx = bart('fakeksp',bart('fftmod -i 6',Imdx),ones(size(params.data)),params.sensemaps);     %multi-coil kspace of dX
FTXFMtdx = bart('fftmod -i 6',FTXFMtdx);
FTXFMtdx=FTXFMtdx.*params.mask;
if params.TVWeight
%     DXFMtx = params.TV*(params.XFM'*x);     %TV OF X
%     DXFMtdx = params.TV*(params.XFM'*dx);   %TV OF DX
else
    DXFMtx = 0;
    DXFMtdx = 0;
end
end


function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params);
%calculated the objective function
p = params.pNorm; 

obj = (FTXFMtx + t*FTXFMtdx) - params.data; %MODIFIED TO ADD PARAM.V
% obj = obj(:)'*(repmat(params.V(:),[size(obj,4) 1]).*obj(:)); %IN THE END THIS IS JUST A VALUE AND NOT A VECTOR= l2-norm
obj = obj(:)'*obj(:); %IN THE END THIS IS JUST A VALUE AND NOT A VECTOR= l2-norm

if params.TVWeight
%     w = DXFMtx(:) + t*DXFMtdx(:);
%     TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end


TV = sum(TV.*params.TVWeight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));
res = obj + (TV) + (XFM) ;
disp(['res: ',num2str(res),' |XFM: ',num2str(XFM),' |obj: ',num2str(obj),' |TV: ',num2str(TV)])

end


function grad = wGradient(x,params)
    
gradXFM = 0;
gradTV = 0;
gradObj = gOBJ(x,params);

if params.xfmWeight
gradXFM = gXFM(x,params);
end

if params.TVWeight
% gradTV = gTV(x,params);
end

grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV);
end %shouldnt change

function grad = gXFM(x,params)
% compute gradient of the L1 transform operator
p = params.pNorm;
grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);
end %shouldnt change

function gradObj = gOBJ(x,params); %gradient of l2-norm

Im=bart('cdf97 -i 6',x);
ksp=bart('fakeksp',bart('fftmod -i 6',Im),params.data,params.sensemaps);        %multi-coil kspace of X
ksp=bart('fftmod -i 6',ksp);
ksp=ksp.*params.mask;
Verror= permute(repmat(params.V,[1 1 size(ksp,4)]),[4 1 2 3]).*(ksp-params.data); %error times variance matrix 
image=bart('fftmod 6',bart('fft 6',bart('fftmod 6',Verror)));
imagec= bart('fmac -C -s8',image, params.sensemaps);
gradObj = bart('cdf97 6',imagec);
gradObj = 2*gradObj; 

end
