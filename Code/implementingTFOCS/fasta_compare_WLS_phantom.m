

Kphantom=fftshift(fft2(phantom(256,256)));
Kphantom=permute(Kphantom,[3 1 2]);


mask=rand(size(Kphantom)); 
mask(1,128-10:128+10,128-10:128+10)=0;
Kphantom((mask>0.5))=0;
[E4,W,Filt,MNSA,mask,pdf,linear_recon_s]=createSPOTops(Kphantom,1);
%%
n1=256; n2=256;
mat = @(x) reshape(x,n1,n2,ncoils);             % function that reshapes vector to matrix 
mat2d = @(x) reshape(x,n1,n2*ncoils);           % functions that reshapes in 2d matrix for visualization purposes 
matcc =@(x) reshape(x,n1,n2);                   % function that reshapes vector to matrix 


KuN=Kphantom(mask); 

MNSAvector=1+floor(10*rand([length(KuN),1]));
MNSAvector=ceil(pdf(mask).*10);
noise=randn([length(KuN),1]).*sqrt(MNSAvector)*1e1; 

KuN=KuN+noise(:);
if false
    MNSA=ones(size(MNSA));
end

mu=1;

disp('FASTA reconstruction')

opts = [];
opts.recordObjective = true;                    % Record the objective function so we can plot it
opts.verbose = true;
opts.stringHeader='    ';                       % Append a tab to all text output from FISTA.  This option makes formatting look a bit nicer. 
opts.accelerate = true;
opts.tol=1e-4; 
% opts.maxIters=50

Wop=opDiag((MNSAvector));

opts.maxIters=100;
muW=mu./(sum(sum(E4*E4'*Wop))./length(KuN)); muW=real(muW);


A=@(x) E4*x;                                    % convert operator to function form 
AT=@(x) E4'*x;


[sol, outs_adapt] = fasta_sparseLeastSquares(A,AT,KuN,mu,ones(size(W*linear_recon_s)), opts);
[solW, outs_adapt] = fasta_sparseweightedLeastSquares(A,AT,Wop,KuN,mu,ones(size(W*linear_recon_s)), opts);

figure(8); imshow(abs([matcc(W'*sol)./max(W'*sol),matcc(W'*solW)./max(W'*solW)]),[])

