function Chrom =MkChrom(nPop,Xres,Zres,acc,Center)
nGenes=ceil(sqrt(ceil(Xres/2)^2+ceil(Xres/2)^2)+1);
Chrom=zeros(nPop,nGenes);




for n=1:nPop
    i=0;
n
while i==0
%     Chrom(n,:)=rand(nGenes,1).^(acc-1)';
        Chrom(n,:)=0.02+rand(nGenes,1).^(acc+2)';

    Chrom(n,1:Center)=ones(1,Center);

    i=CheckAcc(squeeze(Chrom(n,:)),acc);
    end
end