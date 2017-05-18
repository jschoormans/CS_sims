function Chrom =MkChromPolyNom(nPop,Xres,Zres,acc)
nGenes=10;
Chrom=zeros(nPop,nGenes);

for n=1:nPop
   Chrom(n,:)=randn(nGenes,1).*2;
end