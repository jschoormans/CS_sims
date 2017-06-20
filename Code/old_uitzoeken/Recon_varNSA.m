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
    
    mm=length(numbers)
    mk=min(all_klines)
    mak=max(all_klines)

    nchans=max(MR.Parameter.Labels.Index.chan);
    
    V=zeros(mak(1)+abs(mk(1))+1,mak(2)+abs(mk(2))+1,nchans);
    for iline=numbers(MR.Parameter.Labels.Index.typ==1);
        
        chan=MR.Parameter.Labels.Index.chan(iline);
        
        xc=all_klines(iline,1)+1+abs(mk(1));
        yx=all_klines(iline,2)+1+abs(mk(2));
        
        nsa=(V(xc,yx,chan));
        V(xc,yx,chan)=V(xc,yx,chan)+1;
         
        nsa_number(iline) = nsa;
    end
    MR.Parameter.Labels.Index.echo=(nsa_number).';
end

end
end