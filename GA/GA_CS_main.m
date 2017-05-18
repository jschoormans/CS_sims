% QUICK GENETIC ALGO: 
% QUESTION: GIVEN A NOISE CORRUPTED 2D IMAGE, HOW TO DISTRIBUTE SAMLPING POINTS S.T. THE HIGHEST PSNR IS REACHED??
% BEGIN WITH POP OF SAMPLING PATTERNS 
% COST FUNCITON IS FUNCTION OF PSNR
% MUTATE, COMBINE AND ITERATE

% cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/GA')
run GA_CS_init.m
run GA_CS_load_graperfruit


nPop=100; 
nGen=200;
res=size(K,1);
nBest=0.1*nPop;
sig=0;
MutProb=0.001;

acc=3;
P.acc=acc;
P.nsamples=floor(1/acc*res*res);
P.res=res;
P.addcentre=1;
P.usePDF=1;
% MAKE POPULATION OF SAMPLING PATTERNS
Chrom=initGACS(nPop,P);

%% loop over generations
for gen=1:nGen
% EVALUATE save('evolvedmask.mat','MNSA','K','ImRef');

[Fitness,CS,MNSA]=EvalCS(Chrom,K,ImRef,sens,sig) % E is the fitness --> HIGHER = BETTER
%SELECTION 
Selection=CS_GA_SelectAdaptive(Fitness)
% COMBINE 
NewChrom=CS_GA_CrossOver(Chrom,Fitness,Selection,nBest);
% MUTATE
NewChrom=CS_GA_Mutate(NewChrom,MutProb,res,nBest);

Best(gen)=max(Fitness);
MeanFitness(gen)=mean(Fitness);
if gen==1
    CSorig=CS;
end
figure(1); subplot(131);hold on; plot(Best,'b*');plot(MeanFitness,'r*'); hold off;
subplot(232);imshow(CS,[]); axis off
subplot(233);imshow(abs(ImRef),[]);axis off
subplot(235);imshow(MNSA,[0 1]);axis off
subplot(236);imshow(CSorig,[]);axis off



pause(1)
Chrom=NewChrom;
end
%% end loop over generations
export_fig -native '1_EVOL.eps'
export_fig -native '1_EVOL.png'

save('evolvedmask.mat','MNSA','K','ImRef');

%% COMPARE WITH NAIVE MASK
PDF=genPDF([res res],4,1/acc,0,0.15,0);
MNSA_naive1=genSampling(PDF,10,3);
MNSA_naive2=squeeze(bart(['poisson -Y',num2str(res),' -Z',num2str(res),'-y',num2str(acc),'  -z',num2str(acc),' -C5']));

Ki0=makeUndersampledK(K,MNSA,sig);
Ki1=makeUndersampledK(K,MNSA_naive1,sig);
Ki2=makeUndersampledK(K,MNSA_naive2,sig);

CS0=fftshift(abs(bart('pics -l1 -r0.02 -i500',Ki0,ones(size(Ki0)))),2);
CS1=fftshift(bart('pics -l1 -r0.02 -i500',Ki1,ones(size(Ki1))),2);
CS2=fftshift(bart('pics -l1 -r0.02 -i500',Ki2,ones(size(Ki2))),2);

E0=psnr(abs(CS0./max(CS0(:))),abs(ImRef./max(ImRef(:))));
E1=psnr(abs(CS1./max(CS1(:))),abs(ImRef./max(ImRef(:))));
E2=psnr(abs(CS2./max(CS2(:))),abs(ImRef./max(ImRef(:))));


fig2=figure(2);

subplot(231); imshow(MNSA); axis off;xlabel('Evolved mask')
subplot(232); imshow(MNSA_naive1); axis off;xlabel('mask from pdf, p=4')
subplot(233); imshow(MNSA_naive2); axis off;xlabel('mask based on poisson disc')

subplot(234); imshow(abs(CS0),[]); axis off; text(5,5,['PSNR= ',num2str(E0)],'Color','w')
subplot(235); imshow(abs(CS1),[]); axis off; text(5,5,['PSNR= ',num2str(E1)],'Color','white')
subplot(236); imshow(abs(CS2),[]); axis off; text(5,5,['PSNR= ',num2str(E2)],'Color','white')
pos=get(fig2,'position'); set(fig2,'position',[pos(1:2)/4 pos(3:4)*2])
export_fig -native '2_EVOL.eps'
export_fig -native '2_EVOL.png'

%% COMPARE WITH OTHER IMAGE
ImRef2=imread('head256.jpg');
ImRef2=double(ImRef2(:,:,1));
K2=fftshift(fftshift(fft2(ImRef2),1),2);
% K2=K2(128-63:128+64,128-63:128+64);
K2=K2(65:190,65:190);

ImRef2=ifft2(K2);


Ki02=makeUndersampledK(K2,MNSA,sig);
Ki12=makeUndersampledK(K2,MNSA_naive1,sig);
Ki22=makeUndersampledK(K2,MNSA_naive2,sig);

CS02=fftshift(fftshift(bart('pics -l1 -r0.01 -i500',Ki02,ones(size(Ki02))),1),2);
CS12=fftshift(fftshift(bart('pics -l1 -r0.01 -i500',Ki12,ones(size(Ki12))),1),2);
CS22=fftshift(fftshift(bart('pics -l1 -r0.01 -i500',Ki22,ones(size(Ki22))),1),2);

E02=psnr(abs(CS02./max(CS02(:))),abs(ImRef2./max(ImRef2(:))));
E12=psnr(abs(CS12./max(CS12(:))),abs(ImRef2./max(ImRef2(:))));
E22=psnr(abs(CS22./max(CS22(:))),abs(ImRef2./max(ImRef2(:))));


fig3=figure(3);

subplot(231); imshow(MNSA); axis off;xlabel('Evolved mask')
subplot(232); imshow(MNSA_naive1); axis off;xlabel('mask from pdf, p=4')
subplot(233); imshow(MNSA_naive2); axis off;xlabel('mask based on poisson disc')

subplot(234); imshow(abs(CS02),[]); axis off; text(5,5,['PSNR= ',num2str(E02)],'Color','w')
subplot(235); imshow(abs(CS12),[]); axis off; text(5,5,['PSNR= ',num2str(E12)],'Color','white')
subplot(236); imshow(abs(CS22),[]); axis off; text(5,5,['PSNR= ',num2str(E22)],'Color','white')
pos=get(fig3,'position'); set(fig3,'position',[pos(1:2)/4 pos(3:4)*2])
export_fig -native '3_EVOL.eps'
export_fig -native '3_EVOL.png'