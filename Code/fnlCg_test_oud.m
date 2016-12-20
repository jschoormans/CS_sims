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




x = x0;


% line search parameters
maxlsiter = params.lineSearchItnlim ;
gradToll = params.gradToll ;
alpha = params.lineSearchAlpha; ,    beta = params.lineSearchBeta;
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
	[FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params);
	f0 = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, 0, params);
	t = t0;
    [f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
	
	lsiter = 0;

	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:)))^2 & (lsiter<maxlsiter) %FIND T SUCH THAT OBJECTIVE FUNCTION IS MINIMZED ALONG dx
		lsiter = lsiter + 1;
		t = t * beta;
		[f1, ERRobj, RMSerr]  =  objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx,x,dx, t, params);
        
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

	% print some numbers	
    if params.display,
%         fprintf(' ite = %d, cost = %f |lsiter= %f|change=%f \n',k,f1,lsiter,change(k+1));
        
%         figure(param.slice+100);
       subplot(1,2,1)
        
        sf1(k+1)=f1; sL2Obj(k+1)=L2Obj; sL1Obj(k+1)=L1Obj;sL1Obj2(k+1)=L1Obj2;
       
        hold on
        plot([0:k],sf1,'k-')
        plot([0:k],sL2Obj,'g-')
        plot([0:k],param.lambda*sL1Obj,'b-')
        plot([0:k],param.lambda2*sL1Obj2,'r-')
        
        hold off
        legend('f1 (sum)','L2','TV','TGV2')
        drawnow;
        s2=subplot(2,2,2);
        imshow(abs(x0(:,:,floor(size(x0,3)/2))),[]);
        title(s2,'zero-filled recon (one frame)')
        s3=subplot(2,2,4);
        imshow(abs(x(:,:,floor(size(x0,3)/2))),[]);
        title(s3,'CS recon (one frame)');
        drawnow;
        
        %         saveas(gcf,['/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab/Results/temp/temp' filename num2str(param.slice) '.jpg']);
        if ~(k >= param.nite)
            clf; end
        
    end
    %conjugate gradient calculation
    
	g1 = wGradient(x,params);
    if strcmp(params.Beta,'FR')==1  
        bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps)
    elseif strcmp(params.Beta,'PR')==1 
        bk = g1(:)'*(g1(:)-g0(:))/(g0(:)'*g0(:)+eps)
    elseif strcmp(params.Beta,'PR_restart')
        bk = g1(:)'*(g1(:)-g0(:))/(g0(:)'*g0(:)+eps);
        bk=max(real(bk),0); %assume bk should be real??
    elseif strcmp(params.Beta,'LS');
        bk = g1(:)'*(g1(:)-g0(:))/(-dx(:)'*g0(:)+eps);
    elseif strcmp(params.Beta,'LS_restart')
        bk = g1(:)'*(g1(:)-g0(:))/(-dx(:)'*g0(:)+eps);
        bk=max(real(bk),0); %assume bk should be real??
    elseif strcmp(params.Beta,'MLS');
        mu=0.25;
        yk=g1(:)-g0(:);
        bk0 =real(g1(:)'*(yk)/(-dx(:)'*g0(:)+eps))
        bkr=real((mu*(yk(:)'*yk(:))/(-dx(:)'*g0(:)+eps)^2)*g1(:)'*dx(:))
        bk=real(bk0-min(bk0,bkr)) %really not sure about the real part....
    elseif strcmp(params.Beta,'HHS'); %hybrid HS
        bFR = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
        bPR = g1(:)'*(g1(:)-g0(:))/(g0(:)'*g0(:)+eps);
        bk=max(0,min(bFR,bPR));

    end
    
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


function [FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx] = preobjective(x, dx, params)

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





function [res, obj, RMS] = objective(FTXFMtx, FTXFMtdx, DXFMtx, DXFMtdx, x,dx,t, params);
%calculated the objective function

p = params.pNorm; 

obj = params.V.*(FTXFMtx + t*FTXFMtdx) - params.V.*params.data; %MODIFIED TO ADD PARAM.V
% obj = obj(:)'*obj(:); %IN THE END THIS IS JUST A VALUE AND NOT A VECTOR= l2-norm

% obj = (FTXFMtx + t*FTXFMtdx) - params.data; %MODIFIED TO ADD PARAM.V
obj = obj(:)'*(params.V(:).*obj(:)); %IN THE END THIS IS JUST A VALUE AND NOT A VECTOR= l2-norm


if params.TVWeight
    w = DXFMtx(:) + t*DXFMtdx(:);
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
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

function grad = wGradient(x,params)

gradXFM = 0;
gradTV = 0;

gradObj = gOBJ(x,params);
if params.xfmWeight
gradXFM = gXFM(x,params);
end
if params.TVWeight
gradTV = gTV(x,params);
end


grad = (gradObj +  params.xfmWeight.*gradXFM + params.TVWeight.*gradTV);



function gradObj = gOBJ(x,params);
% computes the gradient of the data consistency step: should we here weight
% higher-error terms more (or less)???

% 	gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data));
gradObj = params.XFM*(params.FT'*(params.V.*(params.FT*(params.XFM'*x) - params.data))); %test
% gradObj = params.XFM*(params.FT'*(params.FT*(params.XFM'*x) - params.data)); %test

% no V' : V should be a diagonal matrix, but is reshaped ??
% CHECK V AND FFTSHIFT
if params.Debug==1;
    %% DEBUGGING PLOclosTS
    disp('not sure about debugging functions: check before using')
    Errors=(params.FT*(params.XFM'*x) - params.data);
    ErrorsW=params.V.*(params.FT*(params.XFM'*x) - params.data);
   
   figure(1); hold on; plot(real(Errors(128,:)),'r'); plot(real(ErrorsW(128,:)),'k'); hold off;
   figure(3); hold on; imshow(ErrorsW./(Errors+1e-8),[1 20]); hold off;

[muhat] =  mean(real(Errors(Errors(:)~=0))) 
[sigmahat] =  std(real(Errors(Errors(:)~=0)))   

figure(2); hold on;
[N,X]=hist(real(Errors(Errors(:)~=0)),1000);
hist(real(Errors(Errors(:)~=0)),1000)
norm=(1/(sigmahat*sqrt(2*pi)))*exp((-(X-muhat).^2)./(2*sigmahat.^2));
plot(X,sum(N).*norm./sum(norm))
hold off
end
   %%

gradObj = 2*gradObj ; % DO NOT KNOW WHY IT DOES TWO HERE????

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






