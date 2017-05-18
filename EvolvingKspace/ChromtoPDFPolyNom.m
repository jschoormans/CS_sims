function PDF=ChromtoPDFPolyNom(Chrom,Xres,Zres,acc) %makes 2D PDF From Chrom

L=max(Xres,Zres);
R=linspace(0,1,L);
 A=zeros(Xres,Zres); A(floor(Xres/2)+1,floor(Zres/2)+1)=1; 
 D=round(bwdist(A))+1;

for n=1:size(Chrom,1)
        PolyNom(n,:)=ones(1,L).*Chrom(n,1)+Chrom(n,2).*R+Chrom(n,3).*R.^2+Chrom(n,4).*R.^3+Chrom(n,5).*R.^4+Chrom(n,6).*R.^5+Chrom(n,7).*R.^6+Chrom(n,8).*R.^7+Chrom(n,9).*R.^8;
for r=1:Xres
    PDF(n,r,:)=PolyNom(n,D(r,1:Xres));
end

 PDF(n,:,:)=(Xres*Zres).*PDF(n,:,:)./(acc*sum(sum(PDF(n,:,:),2),3));
end