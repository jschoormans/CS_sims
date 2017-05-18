function S=initGACS(nPop,P)


if P.usePDF==0;
for i =1:nPop
S(i,:,:)=int32(ceil(rand(2,P.nsamples)*P.res));
end


if P.addcentre==1;
    for i =1:nPop
    count=1;
    for j1=floor(P.res/2-5):floor(P.res/2+5);
        for j=floor(P.res/2-5):floor(P.res/2+5);
            
            S(i,1,count)=j1;
            S(i,2,count)=j;
            count=count+1;
            end
        end
        
        
    end
end
    
else %start with pdf based chrom
    %%
    for i =1:nPop
        PDF=genPDF([P.res P.res],3+2*rand(),1/P.acc,0,0,0);
        SamplingMatrix=genSampling(PDF,1,1);
        [index]=find(SamplingMatrix(:));
        [xc,yc]=ind2sub(size(SamplingMatrix),index);
        des_length=floor(length(SamplingMatrix(:))*(1/P.acc));

        ra=randperm(des_length); %permute indices
        xc=xc(ra);
        yc=yc(ra);
        
        if length(xc)>des_length
            xc=xc(1:des_length) ;
            yc=yc(1:des_length) ;
            
        elseif length(xc<des_length)
            
        end
        
        S(i,:,:)=int32([xc,yc]');
    end
    %%
    
end

end