function NewChrom=CrossOver(OldChrom,S,Fitness,IndexC,nBest,acc) %input: oldchrom, selecting matrix
%adaptive crossover


k1=0.5;
k3=k1;
pc=k1*(max(Fitness)-Fitness)./(max(Fitness)-mean(Fitness)); %crossover probability


nPop=size(OldChrom,1);

for n=1:nBest
NewChrom(n,:)=squeeze(OldChrom(IndexC(n),:));
end

for n=nBest+1:nPop
    n
    i=0
    while i==0
   R1=rand(size(OldChrom,2),size(OldChrom,3))<pc(S(1,n));
    R2=~R1;

    NewChrom(n,:)=squeeze(R1.*squeeze(OldChrom(S(1,n),:,:))'+R2.*squeeze(OldChrom(S(2,n),:,:))');

    i=CheckAcc(NewChrom(n,:),acc)
    end
end
end