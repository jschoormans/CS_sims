classdef Recon_varNSA < MRecon
properties
% No additional properties needed
end
methods
    function MR = Recon_varNSA( filename )
        % Create an MRecon object MR upon creation of My_Recon
        MR = MR@MRecon(filename);
    end
    % Overload (overwrite) the existing Perform function of MRecon
    function Perform( MR )
MR.Parameter.Parameter2Read.typ = 1;
% Produce k-space Data (using existing MRecon functions)
MR.ReadData;
MR.addAveragestoLabels
MR.DcOffsetCorrection;
MR.PDACorrection;
MR.RandomPhaseCorrection;
MR.MeasPhaseCorrection;
MR.SortData;
MR.GridData;
MR.RingingFilter;
MR.ZeroFill;
MR.K2IM;
MR.EPIPhaseCorrection;
MR.K2IP;
MR.GridderNormalization;
MR.SENSEUnfold;
MR.PartialFourier;
MR.ConcomitantFieldCorrection;
MR.DivideFlowSegments;
MR.CombineCoils;
MR.Average;
MR.GeometryCorrection;
MR.RemoveOversampling;
MR.ZeroFill;
MR.FlowPhaseCorrection;
MR.RotateImage;

    
end

function addAveragestoLabels(MR)
    disp('unsure about different daata for now: whcih aver label should they have? check in other scan!')
    all_klines=double([MR.Parameter.Labels.Index.ky MR.Parameter.Labels.Index.kz]);
    all_klines(MR.Parameter.Labels.Index.typ~=1)=NaN;
    numbers=1:length(MR.Parameter.Labels.Index.typ);
    for iline=numbers(MR.Parameter.Labels.Index.typ==1);
        loc = find(all(all_klines == ones(size(all_klines,1),1) * all_klines(iline,:),2));
        nsa_number(iline) = find(loc == iline);
    end
    nsa_number=ceil(nsa_number./size(MR.Parameter.Labels.CoilNrs,1)); %temp hack 13 is number of channels!!! MAKE SURE THAT THIS IS CORRECT!!!
    nsa_number=nsa_number-1;
    MR.Parameter.Labels.Index.aver=uint16(nsa_number);
end

end
end