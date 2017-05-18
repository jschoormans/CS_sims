function Chrom=Mutation2(Chrom,MutRate,acc)


MutChrom=rand(size(Chrom))
MutProb=rand(size(Chrom))<MutRate;
Chrom(MutProb)=MutChrom(MutProb);
end

