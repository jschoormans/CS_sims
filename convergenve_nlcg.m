load()

%%
P.Prec=1
P.display=1
P.nite=5;
P.Beta='PR_restart'
tic 
[recon_csPR,costPR] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);
toc
%%
P.Beta='HZ'
tic 
[recon_csHZ,costHZ] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);
toc
%%
P.Beta='FR'
tic 
[recon_csFR,costFR] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);
toc
%%
P.Beta='LS_restart'
tic 
[recon_csLS,costLS] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);
toc
%%
P.Beta='MLS'
tic 
[recon_csMLS,costMLS] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);
toc
%%
P.Beta='HHS'
tic 
[recon_csHHS,costHHS] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);
toc

%%
P.nite=15;
P.Prec=1;
[recon_Prec,costPrec] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);

P.Prec=0
[recon_noPrec,costnoPrec] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);

close all
figure(2);hold on;plot(costPrec);
plot(costnoPrec);hold off; legend(' Preconditioning','No Preconditioning')
figure(3); subplot(121); imshow(abs(recon_noPrec(:,:,2)),[]);
subplot(122); imshow(abs(recon_Prec(:,:,2)),[]);
figure(4); hold on; imshow(squeeze(abs(recon_noPrec(100,:,:))),[])
%%
figure(2);hold on;
plot(costPR);
plot(costHZ);
plot(costFR);
plot(costLS);
plot(costMLS);
plot(costHHS);

hold off
% legend('PR','FR (original)','LS')
legend('PR','FR(orig)','LS','MLS','hybrid HS')
%%
mincost=min([costPR,costFR,costLS])
param.gamma*L1Wgrad+param.alpha*L1TVGrad;
figure(3);
semilogy(costPR-mincost);hold on; 
semilogy(costFR-mincost);semilogy(costLS-mincost);
hold off
legend('PR','FR (original)','LS' )
%%
figure(4);
subplot(221)
imshow(abs(recon_csPR(:,:,2)),[]); axis off

subplot(222)
imshow(abs(recon_csFR(:,:,2)),[]);axis off

subplot(223)
imshow(abs(recon_csLS(:,:,2)),[]);axis off

subplot(224)
imshow(abs(recon_csMLS(:,:,2)),[]);axis off
%%
P.secant=0
P.nite=100;
P.Prec=1;
[recon_Prec,costPrec] = CSL1NlCg_experimental(squeeze(R),P,tempy,tempE,tempslice,filename);

