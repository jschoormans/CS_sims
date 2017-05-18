function [CostS,IndexC,Im,Cost]=Costfunction(Chrom,ImRef,P,Xres,Zres)



nPop=size(Chrom,1)
    parfor n=1:nPop %parfor used to work tho...
        n
        Im(n,:,:)=bart('pics -RW:3:0:0.05 -i25 -S',squeeze(Chrom(n,:,:)).*P,ones(Xres,Zres));
%         Cost(n)=1./sum(sum(abs(squeeze((Im(n,:,:)))-abs(ImRef)).^2,2),1);
        Cost(n)=psnr(abs(squeeze(Im(n,:,:))),abs(ImRef));
    end
    if false
    ImRefN=abs(ImRef)./sum(sum(abs(ImRef).^2,1),2); %NORMALIZE ENERGY
    for n=1:nPop
        ImN(n,:,:)=abs(Im(n,:,:))./sum(sum(abs(squeeze(Im(n,:,:))).^2,1),2);%NORMALIZE ENERGY
            Cost(n)=1./sum((ImRefN(205,:)'-squeeze(ImN(n,205,:))).^2);
    end
    elseif false 
           ImRefN=abs(ImRef)./sum(sum(abs(ImRef).^2,1),2); %NORMALIZE ENERGY
    for n=1:nPop
        ImN(n,:,:)=abs(Im(n,:,:))./sum(sum(abs(squeeze(Im(n,:,:))).^2,1),2);%NORMALIZE ENERGY
        Diff=(diff(abs(squeeze((Im(n,:,:)))))-diff(abs(ImRef))).^2;
        Cost(n)=1./(sum(sum(Diff,2),1))     
    end

    elseif false %focus on one region in image(high detail region)
            for n=1:nPop
%             Cost(n)=1./immse(abs(squeeze(Im(n,100:105,52:76))),abs(squeeze(ImRef(100:105,52:76))))
            Cost(n)=immse(abs(squeeze(Im(n,200:210,105:151))),abs(squeeze(ImRef(200:210,105:151)))).^5
            end 
    elseif false
         for n=1:nPop
            Cost(n)=abs(msssim(abs(squeeze(Im(n,:,:))),abs(squeeze(ImRef(:,:)))))
            end 
    end
    
    [CostS,IndexC]=sort(Cost,'descend')
end