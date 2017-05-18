function [x] = nl_conjgrad_fluor(A, b, x,niter,varargin)

n1=varargin{3}
n2=varargin{4};
rr = @(I) reshape(I,[n1,n2])


T= MakeWaveletOp(x,n1,n2)


if nargin>4
    realI=varargin{1};
    realflag=1;
    lambda=varargin{2}
else
    realflag=0;
    lambda=0; disp('no lambda used')
end


figure(100); 
x=zeros(size(A'*b)); 
subplot(221); imshow(rr(x)); drawnow;

grad=-gradient(A,b,x,T,lambda);
s=grad;


t0=1
bb=0.7
for i=1:niter
    t=t0;
    % perform line search
    
    lsiter=0;
    alpha=1;
    f0=Calcobjective(A,b,x,T,lambda);
    f1=Calcobjective(A,b,x+alpha*s,T,lambda);
    
    
    
    
    while ((f1) > (f0 - alpha*t*abs(s(:)'*grad(:))))  % change this 
        f1=Calcobjective(A,b,x+alpha*s,T,lambda);
        alpha=alpha.*t*bb;
        if lsiter>40; disp('line search failed, convergence reached?');
            
            if realflag==1
                mse=sum((x-realI).^2);
            else
                mse=0;
            end
            
            disp(['iter:',num2str(i),' MSE: ',num2str(mse)])
            
            return; end
        
        lsiter=lsiter+1;
    end
    
    if lsiter > 2
		t0 = t0 * bb;
	end 
	
	if lsiter<1
		t0 = t0 / bb;
	end
    disp(['lsiter=',num2str(lsiter),' alpha=',num2str(alpha)])
    % update the position
    xold=x; 
    x=xold+alpha*s;
    
    % calculate steepest direction
    gradold=grad;
    grad=-gradient(A,b,x,T,lambda);
    
    % calculate beta: TO DO 
    beta=(grad'*(grad-gradold))/(gradold'*gradold+eps); %POLAK-RIBIERE 
%     beta=0 ;% Newton??
    
    % update conjugate direction:
    sold=s;
    s=grad+beta*sold;
    
    subplot(221);
    imshow(rr(abs(x)),[]); drawnow;
    if realflag==1
    mse=sum((x-realI).^2);
    else
        mse=0;
    end
    pause(1);
    
    
    disp(['iter:',num2str(i),' MSE: ',num2str(mse)])

end



end

function grad=gradient(A,b,x,T,lambda) % for image domain at least; 

gradl1=gradientl1(A,b,x,T);
grad=2*A'*(A*x-b)+lambda*gradl1; 

end

function objective=Calcobjective(A,b,s,T,lambda)
obj=(A*s-b);
objectivel2=obj(:)'*obj(:);
objectivel1=objectivel1(A,b,s,T);

objective=sum(objectivel2+lambda*sum(objectivel1(:)));
end


function objective=objectivel1(A,b,s,T)
objective=((T*s).*conj(T*s)+eps).^(1/2); 
end

function gradl1=gradientl1(A,b,x,T)
x=T*x; 
gradl1 = x.*(x.*conj(x)+eps).^(-1/2);

end

function W=MakeWaveletOp(x,n1,n2)

W=opWavelet2(n1,n2,'Daubechies',4,4,0);

end
