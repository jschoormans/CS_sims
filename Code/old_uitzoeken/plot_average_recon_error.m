% PLOT AVERAGE RECON ERROR

sigma2=1e-5; %fix
n=1000; %fix
z=randn(1,n);
us=linspace(0,1,100);
k=linspace(0,1,250).*n;%sparsity


for i=1:100
    for j=1:length(k);
        undersampling=us(i);
        sparsity=k(j);
        zu=linspace(0,1,length(z))>(1-undersampling);
        F(i,j)=sparsity*sigma2*log(n)*sum(sqrt(abs(zu.*z).*abs(zu.*z)));
    end
end

%
figure(1); surf(k,us,F)
ylabel('undersampling')
xlabel('sparsity')
zlabel('average recon error')
%% make sigma2 a function of us 
clear F
n=1000; %fix
z=randn(1,n);
us=linspace(0,1,100);
k=linspace(0,1,250).*n;%sparsity


for i=1:100
    for j=1:length(k);


        undersampling=us(i);
        sparsity=k(j);
        sigma2=1e-5/sqrt(1/undersampling+eps);
        
        zu=linspace(0,1,length(z))>(1-undersampling);
        F(i,j)=sparsity*sigma2*log(n)*sum(sqrt(abs(zu.*z).*abs(zu.*z)));
    end
end
%
figure(2); surf(k,us,F)
ylabel('undersampling/averaging')
xlabel('sparsity')

%%
figure(3)
hold on 
plot(us,F(:,10))

plot(us,F(:,30))

plot(us,F(:,50))

plot(us,F(:,100))
plot(us,F(:,200))
hold off
xlabel('undersampling/averaging')
ylabel('average recon error')
legend(' sparse',' less sparse')
%% WHAT IS THE VALUE OF ADDING MEASUREMTENS (REMOVING NOISE??)
clear F
n=100; %fix
z=randn(1,n);
us=linspace(0,1,100);
k=linspace(0,1,25).*n;%sparsity
nmeas=10;

for i=1:nmeas %
    for j=1:length(k);
        undersampling=0.9;
        sparsity=k(j);
        sigma2=1e-5/sqrt(i)
        
        zu=linspace(0,1,length(z))>(1-undersampling);
        F(i,j)=sparsity*sigma2*log(n)*sum(sqrt(abs(zu.*z).*abs(zu.*z)));
    end
end
%
figure(4); plot([2:10],cumsum(-diff(F)))
xlabel('NSA')
ylabel('decrease in average reconstruction error')
title('What is the value of adding measurements? sparsity 0 to 100 pct')


