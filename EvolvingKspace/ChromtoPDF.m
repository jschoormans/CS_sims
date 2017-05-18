function PDF=ChromtoPDF(Chrom,Xres,Zres) %makes 2D PDF From Chrom

 A=zeros(Xres,Zres); A(floor(Xres/2)+1,floor(Zres/2)+1)=1; 
 D=round(bwdist(A))+1;

for n=1:size(Chrom,1)
for r=1:Xres
    PDF(n,r,:)=Chrom(n,D(r,:));
end
end