function NewChrom=CS_GA_CrossOver(S,Fitness,Selection,nBest)
%adaptive crossover


k1=0.5;
k3=k1;
pc=k1*(max(Fitness)-Fitness)./(max(Fitness)-mean(Fitness)); %crossover probability

nPop=size(S,1);

% TO DO: ADD SURVIVING BEST INDS


[~,i]=sort(Fitness,'descend');
for n=1:nBest
   NewChrom(n,:,:) =S(i(n),:,:);
end
    
for n=nBest+1:nPop
    R1=repmat(rand(1,size(S,3))>pc(Selection(1,n)),[2,1]);
    R2=~R1;
    R1=int32(R1); R2=int32(R2);
    NewChrom(n,:,:)=squeeze(S(Selection(1,n),:,:)).*R1+squeeze(S(Selection(2,n),:,:)).*R2;
end
end