mu = .1;
lambda = mu/5;
delta=mu;
gamma = mu/1000;
values=[mu,mu/5,mu/10,mu/50,mu/100];
kkx=250;
Rtemp=squeeze(R(kkx,:,:));
Datatemp=squeeze(Data(kkx,:,:));
parfor idelta=1:5
idelta
for ilambda=1:5
delta=values(idelta); lambda=values(ilambda);
u{idelta,ilambda}= mricswaveletTV(Rtemp,Datatemp, mu, lambda, gamma,delta, nInner, nBreg);
end
end

%%
for idelta=1:5
idelta
for ilambda=1:5
W(idelta,ilambda)=sum(sum(abs(fwt2(u{idelta,ilambda},'db4',4))))
f(idelta,ilambda)=sum(sum((fft2(u{idelta,ilambda})-Datatemp).^2))
end
end
%%
figure(1)
for idelta=1:5
for ilambda=1:5
subplot(5,5,(idelta-1)*5+ilambda)
imshow(abs(u{idelta,ilambda}),[])
title(strcat('\lambda = ',num2str(values(ilambda)),' & \delta = ',num2str(values(idelta))))
end
end

%%
figure(2)

surf(W)
xlabel('delta')
ylabel('lambda values')

%%
figure(3)
surf(abs(f))
xlabel('delta')
ylabel('lambda values')