function M=MutationAdaptive(Chrom,MutRate,acc,adaptiveacc)

if adaptiveacc==1
   acc=sum(Chrom(:))./(size(Chrom,1)*size(Chrom,2)) 
end




if true
Mutplus=rand(size(Chrom))<MutRate;  %mutates zeroes into ones
Mutmin=rand(size(Chrom))<MutRate*acc;  %mutates ones into zeroes 

% do M-: mutate ones to zeroes; zeroes stay zeros;
M2=Chrom>Mutmin;

% do M+; mutates zeroes to ones; ones stap ones
M3=M2+(Mutplus>M2);
M=M3;

else % mutate points with pdf
%possibilty?/////////////////
end

