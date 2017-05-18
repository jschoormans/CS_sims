function S=SelectAdaptive(Fitness,Index,nBest) %input: oldchrom, selecting matrix

nPop=size(Index,2);

    S(1,:)=sus(Fitness',nPop)
    S(2,:)=sus(Fitness',nPop)
    
end