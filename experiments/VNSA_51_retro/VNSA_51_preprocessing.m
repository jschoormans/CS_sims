% VNSA_50 preprocessing: 
% GOAL: CONVERT .RAW in .mat KSPACES

folder='/home/jschoormans/lood_storage/divi/Ima/parrec/Jasper/VNSA/VNSA_51/';
cd(folder);

F=dir('*.raw')

for ii=[4,5]
    cd(folder);
    file=F(ii).name
    
    MR=MRecon([folder,file]);
    MR.Parameter.Parameter2Read.typ=1;
    MR.Parameter.Recon.ImmediateAveraging='no'
    MR.ReadData;
    MR.RandomPhaseCorrection;
    MR.PDACorrection;
    MR.DcOffsetCorrection;
    MR.MeasPhaseCorrection;
    MR.SortData;
    % MR.Average;
    MR.RemoveOversampling;
    MR.RingingFilter;
    MR.ZeroFill;
    % MR.K2I;
    % MR.ShowData
    K=MR.Data;
    K=squeeze(K);
    K=K./max(K(:)); %normalize kspace
    
    cd('/scratch/jschoormans/VNSA_51')
    filename = strtok(file,'.')
    save(filename,'K','-v7.3');
    
end