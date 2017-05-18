function NewChrom=CS_GA_Mutate(S,MutProb,res,nBest)
%adaptive mutation
nPop=size(S,1);
P.nsamples=size(S,3);
P.res=res;
P.addcentre=0;
P.usePDF=0;
SNew=initGACS(nPop,P);
for n=1:nBest
       NewChrom(n,:,:)=S(n,:,:);
end
for n=nBest+1:nPop
    if true %rando mutataion
    R1=rand(1,2,P.nsamples)<MutProb;
    R2=~R1;
    NewChrom(n,:,:)=int32(R2).*S(n,:,:)+int32(R1).*SNew(n,:,:);
    end
    if true %mutation by moving around 
        R1=rand(1,2,P.nsamples)<MutProb;
        R2=~R1;
        O=(-1.*int32(ones(1,2,P.nsamples))+2*int32(rand(1,2,P.nsamples)>0.5));
        NewChrom(n,:,:)=int32(R2).*S(n,:,:)+int32(R1).*(S(n,:,:)+O);
    end
end

if true %check wrong integers
    NewChrom(NewChrom>P.res)=1;
    NewChrom(NewChrom==0)=P.res;
end
end