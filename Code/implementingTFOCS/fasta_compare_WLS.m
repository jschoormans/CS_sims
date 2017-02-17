

KuN=Kphantom(:); 

if false
    MNSA=ones(size(MNSA));
end
mu=0.001;

Wop=opDiag(repmat(MNSA(mask).',[ncoils 1]));
opts.maxIters=100;

muW=mu./(sum(sum(E4*E4'*Wop))./length(Ku)); muW=real(muW);


[sol, outs_adapt] = fasta_sparseLeastSquares(A,AT,KuN,mu,ones(size(W*linear_recon_s)), opts);
[solW, outs_adapt] = fasta_sparseweightedLeastSquares(A,AT,Wop,KuN,muW,ones(size(W*linear_recon_s)), opts);

figure(8); imshow(abs([matcc(W'*sol)./max(W'*sol),matcc(W'*solW)./max(W'*solW), MNSA./max(MNSA(:))]),[])

