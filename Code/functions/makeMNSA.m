function MNSA= makeMNSA(P)

MNSA =zeros(size(P.pdf));
epsilon=0;
% make weighting scheme
if P.jjj==1 %SHOULD BE CHANGED
    MNSA=ceil(1./P.pdf);
elseif P.jjj==2 %SHOULD BE CHANGED
    MNSA=ceil(P.pdf*P.acc);
    
elseif P.jjj==3
    MNSA=P.acc*ones(size(P.pdf));
elseif P.jjj==4
    MNSA=ones(size(P.pdf));
elseif P.jjj==5 %USEDYNS FROM HERE
    MNSA=P.acc*ones(size(P.pdf))*P.usedyns;
elseif P.jjj==6 % MORE OUTSIDE
    MNSA=ceil(1./P.pdf)*P.usedyns;
elseif P.jjj==7 % MORE IN CENTER
    MNSAref=P.acc*ones(size(P.pdf))*P.usedyns;
    refsum=sum(sum(MNSAref.*P.M))
    epsilon=refsum./round(sum(sum(P.pdf.*P.acc.*P.M)));    
    MNSA=round(epsilon*(P.pdf*P.acc));

elseif P.jjj==8 %random weighting of points
    MNSAref=P.acc*ones(size(P.pdf))*P.usedyns;
    refsum=sum(sum(MNSAref.*P.M));
    epsilon=1
        MNSA=floor(epsilon*P.acc*rand(size(P.pdf))*P.usedyns*2)+1;

    while abs(sum(sum(MNSA.*P.M))-refsum) <2e1
    MNSA=floor(epsilon*P.acc*rand(size(P.pdf))*P.usedyns*2)+1;
    sum(sum(MNSA.*P.M))
    epsilon=refsum/sum(sum(MNSA.*P.M))
    end
elseif P.jjj==9 % MORE IN CENTER
    MNSAref=P.acc*ones(size(P.pdf))*P.usedyns;
    refsum=sum(sum(MNSAref.*P.M))
    epsilon=refsum./round(sum(sum((P.pdf.^(0.8)).*P.acc.*P.M)));    
    MNSA=round(epsilon*((P.pdf.^(0.8)).*P.acc));
end
    
end

