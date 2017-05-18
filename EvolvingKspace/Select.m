function S=Select(IndexC,nBest) %input: oldchrom, selecting matrix

nPop=size(IndexC,2);

for n=1:nBest %the best ones should be always selected
S(1,n)=IndexC(n);
S(2,n)=IndexC(n);

end
for n=nBest+1:nPop
    
    %cumulative distirubtion function: one is most likely to be chosen;
    r1=SelectRandom(nPop)
    r2=SelectRandom(nPop)
    
    S(1,n)=IndexC(r1)
    S(2,n)=IndexC(r2)
end

end