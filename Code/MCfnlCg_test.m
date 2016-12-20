function x = MCfnlCg_test(x0,params)
%-----------------------------------------------------------------------
%
% res = fnlCg(x0,params)
%
% implementation of a L1 penalized non linear conjugate gradient reconstruction
%
% The function solves the following problem:
%
% given k-space measurments y, and a fourier operator F the function 
% finds the image x that minimizes:
%
% Phi(x) = ||F* W' *x - y||^2 + lambda1*|x|_1 + lambda2*TV(W'*x) 
%
%
% the optimization method used is non linear conjugate gradient with fast&cheap backtracking
% line-search.
% 
% (c) Michael Lustig 2007
%-------------------------------------------------------------------------

% TO DO: 
% ADD MULTICOIL OPERATOR 
% ADD BETTER (FASTER?) CG ALGO?


x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha;
beta = params.lineSearchBeta;
t0 = params.lineSearchT0;
k = 0;
t = 1;

% copmute g0  = grad(Phi(x))

g0 = wGradient(x,params); %x is wavelet image; 

dx = -g0;


% iterations
while(1)

% backtracking line-search

	% pre-calculate values, such that it would be cheap to compute the objective
	% many times for efficient line-search
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,DXFM2tx, DXFM2tdx] = preobjective(x, dx, params);
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,DXFM2tx, DXFM2tdx,x,dx, 0, params);
	t = t0;
    [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,DXFM2tx, DXFM2tdx,x,dx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter) %FIND T SUCH THAT OBJECTIVE FUNCTION IS MINIMZED ALONG dx
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,DXFM2tx, DXFM2tdx,x,dx, t, params);
	end

	if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	% control the number of line searches by adapting the initial step search
	if lsiter > 2
		t0 = t0 * beta;
	end 
	
	if lsiter<1
		t0 = t0 / beta;
	end

	x = (x + t*dx);

	%--------- uncomment for debug purposes ------------------------	
	disp(sprintf('%d   , obj: %f, RMS: %f, L-S: %d, ErrObj: %f', k,f1,RMSerr,lsiter,ERRobj));
    if params.display==1
        figure(100);  
        subplot(121)
        sf1(k+1)=abs(f1);
        plot([0:k],sf1,'k-');
        subplot(222);
        imshow(params.V(:,:,1),[]);axis off;
        subplot(224);
        imshow(abs(params.XFM'*x),[]);axis off;

        drawnow; 
    end
	%---------------------------------------------------------------
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx =  - g1 + bk* dx;
	k = k + 1;
	
	%TODO: need to "think" of a "better" stopping criteria ;-)
	if (k > params.Itnlim) | (norm(dx(:)) < gradToll) 
		break;
	end

end


return;


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,DXFM2tx, DXFM2tdx] = preobjective(x, dx, params)

% precalculates transforms to make line search cheap

FTXFMtx = params.FT*(params.XFM'*x);
FTXFMtdx = params.FT*(params.XFM'*dx);

if params.TVWeight
    DXFMtx = params.TV*(params.XFM'*x);
    DXFMtdx = params.TV*(params.XFM'*dx);
else
    DXFMtx = 0;
    DXFMtdx = 0;
end
if params.TV2Weight
    DXFM2tx = params.TV2*(params.XFM'*x);
    DXFM2tdx = params.TV2*(params.XFM'*dx);
else
    DXFM2tx = 0;
    DXFM2tdx = 0;
end



function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,DXFM2tx, DXFM2tdx, x,dx,t, params);
%calculated the objective function

p = params.pNorm; 

% obj = params.V.*(FTXFMtx + t*FTXFMtdx) - params.V.*params.data; %MODIFIED TO ADD PARAM.V
% obj = obj(:)'*obj(:); %IN THE END THIS IS JUST A VALUE AND NOT A VECTOR= l2-norm

obj = (FTXFMtx + t*FTXFMtdx) - params.data; %MODIFIED TO ADD PARAM.V
obj = obj(:)'*(params.V(:).*obj(:)); %IN THE END THIS IS JUST A VALUE AND NOT A VECTOR= l2-norm


if params.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.TV2Weight
    w = DXFM2tx(:) + t*DXFM2tdx(:);
    TV2 = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV2 = 0;
end

if params.xfmWeight
   w = x(:) + t*dx(:); 
   XFM = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    XFM=0;
end



TV = sum(TV.*params.TVWeight(:));
TV2 = sum(TV2.*params.TV2Weight(:));
XFM = sum(XFM.*params.xfmWeight(:));
RMS = sqrt(obj/sum(abs(params.data(:))>0));

res = obj + (TV) + (XFM) ;

function grad = wGradient(x,params)

gradXFM = 0;
gradTV = 0;
gradTV2 = 0;

gradObj = gOBJ(x,params);
if params.xfmWeight
gradXFM = gXFM(x,params);
end
if params.TVWeight
gradTV = gTV(x,params);
end
if params.TV2Weight
gradTV2 = gTV2(x,params);
end

grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV+params.TV2Weight.*gradTV2);

function gradObj = gOBJ(x,params);
% computes the gradient of the data consistency step: should we here weight
% higher-error terms more (or less)???

% 	gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data));
gradObj = params.XFM*(params.FT'*(params.V.*(params.FT*(params.XFM'*x) - params.data))); %test

% no V' : V should be a diagonal matrix, but is reshaped ??
% CHECK V AND FFTSHIFT

   %%

gradObj = 2*gradObj ; 

function grad = gXFM(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

grad = p*x.*(x.*conj(x)+params.l1Smooth).^(p/2-1);

function grad = gTV(x,params)
% compute gradient of TV operator
p = params.pNorm;
Dx = params.TV*(params.XFM'*x);
G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.XFM*(params.TV'*G);


function grad = gTV2(x,params)
% compute gradient of TV operator
p = params.pNorm;
Dx = params.TV2*(params.XFM'*x);
G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.XFM*(params.TV2'*G);


