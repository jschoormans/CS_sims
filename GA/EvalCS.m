function [E,bestCS,bestMNSA]=EvalCS(S,K,ImRef,sens,sig)
nPop=size(S,1);
N=randn(size(K)).*sig; %noise contribution to k-space;
XFM = Wavelet('Daubechies',4,4);	% Wavelet\

poolobj=gcp;
addAttachedFiles(poolobj,{'/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/sparseMRI_v0.2/@p2DFT','/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/sparseMRI_v0.2/@Wavelet','/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab/Misc/@Wavelet/private'})
parfor jj=1:(nPop)
Si=S(jj,:,:);
[E(jj),CS{jj},MNSA{jj}]=GA_CS_CS(Si,K,ImRef,sens,N,XFM);
jj;
end


[~,ij]=sort(E,'descend');
bestMNSA=MNSA{ij(1)};
bestCS=CS{ij(1)};
end