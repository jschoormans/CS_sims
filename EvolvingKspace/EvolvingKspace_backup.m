%EVOLVING KSPACE PATTERNS
clear all; close all; clc;
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/EvolvingKspace')

Xres=256;
Zres=256;
acc=8;

P=bart('phantom -x256 -k');
M=bart('poisson -Y256 -Z256');
ImRef=bart('pics -RW:3:0:0.1 -i25 -S',P,ones(Xres,Zres));

%%
PDF=genPDF([Xres Zres],4,1/acc);
S=genSampling(PDF,10,1);

CS=bart('pics -RW:3:0:0.05 -i25 -S',S.*P,ones(Xres,Zres));
figure(10); subplot(131); imshow(abs(ImRef),[]); title('ref im'); subplot(132); imshow(abs(CS),[]); title('CS with typical pdf')
sum(sum((squeeze(abs(CS))-abs(ImRef)).^2,2),1)

%% PARAMS
clearvars Im Cost BestCost MeanCost Chrom
nGen=5e2; nPop=2e2;
nBest=10;
MutRate=1e-3;
tol=5e-2
p=2;      %determines cost function: higher p means best chroms will be selected more often

%% MAKE CHROM (INITAL POPULATION)
Chrom=zeros(nGen,nPop,Xres,Zres);
for n=1:nPop
    C=0
    PDF=genPDF_nocheck([Xres Zres],3+(10*rand).^2,1/acc,1,C);
    Chrom(1,n,:,:)=genSampling(PDF,10,1);
end
% Chrom(1,:,:,:)=rand(nPop,Xres,Zres)>(1-1/acc);

%% GENETIC LOOP
for gen=1:nGen
   gen
   [CostS,IndexC,Im,Cost]=Costfunction(Chrom(gen,:,:,:),ImRef,P,Xres,Zres);

    BestCost(gen)=CostS(1);
    MeanCost(gen)=mean(CostS);
    
    %show best image of generation
    
    figure(1)
    subplot(131)
    imshow(squeeze(abs(Im(IndexC(1),:,:))),[])
    subplot(132)
    imshow(squeeze(abs(Chrom(gen,IndexC(1),:,:))),[])
    subplot(133)
    imshow(squeeze(mean(abs(Chrom(gen,:,:,:)),2)),[0 1])
    drawnow;
    
    figure(2)
    hold on 
    plot(BestCost,'b.-'); drawnow;
    plot(MeanCost,'r.-'); drawnow;
    hold off
    
    
    %make new generation of masks
    S=SelectAdaptive(CostS,IndexC,nBest)
    Chrom(gen+1,:,:,:)=CrossOver(squeeze(Chrom(gen,:,:,:)),S,Cost,IndexC,nBest);
    
    
    
    for n=nBest+1:nPop
        Chrom(gen+1,n,:,:)=Mutation(squeeze(Chrom(gen+1,n,:,:)),MutRate,acc,1);
        
        if sum(sum(Chrom(gen+1,n,:,:),3),4)>(Zres*Xres*(1/acc)*(1+tol));
            %remove a few points randomly
            while sum(sum(Chrom(gen+1,n,:,:),3),4)>(Zres*Xres*(1/acc)*(1+tol));
                Chrom(gen+1,n,:,:)=(squeeze(Chrom(gen+1,n,:,:)) & rand(Zres,Xres)>0.01);
                disp('removing points...');
                sum(sum(Chrom(gen+1,n,:,:),3),4);
            end
            
        end
    end

end