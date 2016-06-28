classdef CS_Recon < MRecon
    properties
        %none 
    end
    
    methods
        function MR = CS_Recon( filename )
            MR=MR@MRecon(filename);
        end
        % Overload (overwrite) the existing Perform function of MRecon    
        function Perform1( MR )            
            %Reconstruct only standard (imaging) data
            MR.Parameter.Parameter2Read.typ = 1;                        
            % Produce k-space Data (using MRecon functions)
            disp('Reading data...')
            MR.ReadData;
            disp('Corrections')
            MR.DcOffsetCorrection;
            MR.PDACorrection;
            MR.RandomPhaseCorrection;
            MR.MeasPhaseCorrection;
            disp('Sorting data...')
        end
        function Sort(MR)
            MR.SortData;
            disp('Perform part 1 finished')
        end
        function Perform2(MR)
            %disp('Ringing Filter')
            %MR.RingingFilter;
            MR.ZeroFill
            disp('klopt dit niet eerst k2ip??') 
            MR.K2IM %if reconstructing slice by slice (yz) first iFFT in slice-direction
            MR.EPIPhaseCorrection; %EPI correction for FOV/2 Ghost from eddy current effects
            MR.K2IP;
            MR.GridderNormalization;
            disp('Reconstruction done, last few things now...')
            MR.ConcomitantFieldCorrection;
            MR.CombineCoils;
            MR.Average;
            MR.GeometryCorrection;
            MR.RemoveOversampling;
            MR.ZeroFill;
            MR.RotateImage;
            disp('Finished')
        end

    end
    
    % These functions are Hidden to the user
    methods (Static, Hidden)

    end
end