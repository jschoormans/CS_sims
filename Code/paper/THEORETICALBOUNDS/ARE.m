%optimize AvReconError

N=512;
r=linspace(0,1,100);

k=abs(1-r).*N; %sparsity
sigma1=1e-4
undersampling=0.2

z=randn(1,n);
zu=linspace(0,1,length(z))>(1-undersampling);



%%
N(1,:)=2.*ones(size(r));
N(2,:)=ones(size(r))+2.*abs(1-r);
N(3,:)=ones(size(r))+2.*abs(r);

sum(N.')
for i=1:3
    
ARE(i)=sum(k.*(sigma2./sqrt(N(i,:)+eps)))*log(n)*sum(sqrt(abs(zu.*z).*abs(zu.*z)))
end