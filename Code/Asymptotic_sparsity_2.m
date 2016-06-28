%Asymptotic sparsity pt II

% according to Hansen: (slide 49)
%m >~ log(N) * s
%where s is the sparsity (Number of relevant coeffs)
%N is the total number of voxels/coeffs
% m: the number of measurements


run Asymptotic_sparsity.m
%%

epsilon=0.8;
for ii=2:nmax
    sparsity(ii)=(find(F{ii}.C>epsilon,1,'first')./length(F{ii}.C)) ;
    s(ii)=sparsity(ii)*length(F{ii}.C);
    % PDF(2^(ii-1):2^ii)=sparsity(ii);
end

for ii=1:nmax
N(ii)=2^(ii);
end

sc=cumsum(s) %cumulative sum;
%%
figure(10);

subplot(221)
hold on
plot(N,sc,'k')
plot(N,N.^2,'r--')
hold off
title('sparsity')
xlabel('N')
ylabel('s')

subplot(222)
hold on
plot(N,log(N),'k')
hold off
title('log(N)')
xlabel('N')
ylabel('log(N)')



subplot(223)
hold on
plot(N,log(N).*sc,'k')
plot(N,N.^2,'r--')
hold off
title('m >~ log(N)*sc')
ylabel('m')
xlabel('N')

subplot(224)
hold on
plot(N,log(N).*sc./(N.^2),'k')
hold off
title('m/N^2')
xlabel('N')
ylabel('m')

%%%%%%%%%%%%%%%%%%%%%
%% NUMBER OF SAMPLES IN EACH LEVEL
epsilon=0.1;
log(epsilon)