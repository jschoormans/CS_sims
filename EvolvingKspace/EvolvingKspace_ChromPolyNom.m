%EVOLVING KSPACE PATTERNS
clear all; close all; clc;
cd('/home/jschoormans/lood_storage/divi/Projects/cosart/CS_simulations/EvolvingKspace')
addpath(genpath('/home/jschoormans/lood_storage/divi/Projects/cosart/Matlab_Collection/exportfig'))

C=clock;
experiment=[date,'-',num2str(C(4)),'-',num2str(C(5))];
direxp=[cd,'/Experiments/',experiment]
dir=mkdir(direxp)

Xres=256;
Zres=256;
acc=6;

% KSPACE
P=bart(['phantom -x',num2str(Xres),' -k']); %k-space
Im=imread('brain.png');

% P=ifftshift(ifftshift((fft2(double(squeeze(Im(1:512,1:512,1))))),2),1);
Xres=size(P,1);
Zres=size(P,2);




% M=bart(['poisson -Y',num2str(Xres),' -Z',num2str(Zres)]);
ImRef=bart('pics -RW:3:0:0.1 -i25 -S',P,ones(Xres,Zres));

%
PDF=genPDF([Xres Zres],4,1/acc);
S=genSampling(PDF,10,1);

CS=bart('pics -RW:3:0:0.05 -i25 -S',S.*P,ones(Xres,Zres));
figure(10); subplot(131); imshow(abs(ImRef),[]); title('ref im'); subplot(132); imshow(abs(CS),[]); title('CS with typical pdf')
cd(direxp);
export_fig '10.tiff' -native

%% PARAMS
clearvars Im Cost BestCost MeanCost Chrom
nGen=5e2; nPop=4e1;
nBest=10;
MutRate=1e-2;
tol=5e-2
p=2;      %determines cost function: higher p means best chroms will be selected more often

%% MAKE CHROM (INITAL POPULATION)
% Chrom=MkChromPolyNom(nPop,Xres,Zres,acc);
Chrom=MkChrom(nPop,Xres,Zres,acc,10);

% Chrom(1,:,:,:)=rand(nPop,Xres,Zres)>(1-1/acc);

% GENETIC LOOP
for gen=1:nGen
   gen
   
   PDF=ChromtoPDF(Chrom,Xres,Zres);
%       PDF=ChromtoPDFPolyNom(Chrom,Xres,Zres,acc);

   for n=1:nPop
       n
   SamplingPattern(n,:,:)=genSampling(squeeze(PDF(n,:,:)),10,30);
   end
   [CostS,IndexC,Im,Cost]=Costfunction(SamplingPattern,ImRef,P,Xres,Zres);

    BestCost(gen)=CostS(1);
    MeanCost(gen)=mean(CostS);
    
    %show best image of generation
    
    figure(1)
    subplot(131)
    imshow(squeeze(abs(Im(IndexC(1),:,:))),[])
    subplot(232)
    imshow(squeeze(abs(SamplingPattern(IndexC(1),:,:))),[])
    subplot(133)
    imshow(squeeze(mean(abs(SamplingPattern(:,:,:)),1)),[0 1])
    
        subplot(235)
plot(Chrom(IndexC(1),:))
    drawnow;
export_fig '1.tiff' -native

    
    figure(2)
    hold on 
    plot(BestCost,'b.-'); drawnow;
    plot(MeanCost,'r.-'); drawnow;
    hold off
export_fig '2.tiff' -native

    
    %make new generation of masks
    S=SelectAdaptive(CostS,IndexC,nBest)
    NewChrom(:,:)=CrossOver(squeeze(Chrom(:,:)),S,Cost,IndexC,nBest,acc);
    
    NewChrom(nBest+1,:)=Mutation2(NewChrom(nBest+1,:),MutRate,acc);

    check=CheckAcc(Chrom,acc);
    for n=nBest+1:nPop
        while check(n)~=1
            Chrom(n,:)=Chrom(n,:).*0.98;
             check=CheckAcc(Chrom,acc);
        end
    end
    
    Chrom=NewChrom;
end