function S=CS_GA_SelectAdaptive(Fitness) %input: oldchrom, selecting matrix

nPop=size(Fitness,2);

S(1,:)=sus(Fitness',nPop)
S(2,:)=sus(Fitness',nPop)

end