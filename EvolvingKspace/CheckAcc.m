function check=CheckAcc(Chrom,acc);
%checks chromosomes, if they are still good (not too many sampling points);

m=zeros(size(Chrom));
m(:,floor(length(Chrom)/2)+1)=ones(size(Chrom,1),1);
r=bwdist(m);
o=2*pi.*r;

check=sum(Chrom.*o,2).*acc<(sum(o,2));