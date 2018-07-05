
addpath('.\..\..\')
addpath('.\..\..\Tools')
addpath('.\..\..\Methods')

% close all
clear all
%Ask for the output directory and the folders with the input data
DATADIRS = uipickfiles();
for dd = 1:length(DATADIRS)
    DATADIR = [DATADIRS{dd},'\'];
    WaveDirections = {'Upstream','Downstream'};
    
    for ww = 1:length(WaveDirections) 
        WaveDirection = WaveDirections{ww};
        if strcmp(WaveDirection,'Upstream')
            load([DATADIR,'ComplexPressures_Source1.mat'])
            load([DATADIR,'MeasurementInfo.mat'])
            SortedFrequencyIndex = Measurement.fIndex.Rep001(:,1);
        elseif  strcmp(WaveDirection,'Downstream')
            load([DATADIR,'ComplexPressures_Source2.mat'])
            load([DATADIR,'MeasurementInfo.mat'])
            SortedFrequencyIndex = Measurement.fIndex.Rep001(:,2);
        else
            error('WaveDirection incorrectly defined')
        end
        
        Z_cal  = load('.\Calibration\Cal_2018_04_10.mat');
        C_fitted = zeros(size(Z_cal.C,2),length(MeasData.Ref_Freq));
        for ii=1:size(Z_cal.C,2)    
            C_fitted(ii,:) = interp1(Z_cal.f_cal, real(Z_cal.C(:,ii)),MeasData.Ref_Freq,'linear','extrap') + 1i*interp1(Z_cal.f_cal, imag(Z_cal.C(:,ii)),MeasData.Ref_Freq,'linear','extrap');
        end        
        MeasData.Z = MeasData.Z./C_fitted; 
        Epsilon = 65e-2;
        X = 0:0.055:9*0.055;

        MicPositions = X.';
        MicEqPositions = X.';  
        P = MeasData.Z(4:end-3, SortedFrequencyIndex );
        %Calibrate the pressures
        
        if strcmp(WaveDirection,'Upstream')
            FlowVelocity = Measurement.U*Measurement.WaveCal.Velocity_Corr.Port2;
        elseif  strcmp(WaveDirection,'Downstream')
            FlowVelocity = Measurement.U*Measurement.WaveCal.Velocity_Corr.Port2;
        else
            error('WaveDirection incorrectly defined')
        end
        Temperature = Measurement.t;
        Height = 0.025;
        Prony = PronyImpedance( Measurement.f.',P,MicPositions,MicEqPositions,5,Epsilon,...
                                FlowVelocity, Temperature, Height);


        Prony.CalculateImpedance;
        Prony.PlotImpedance;
        %Stop executing and see if the data is correct, use dbcont to save
        %the data and go to the next datafile
        keyboard
        save([DATADIR,WaveDirection,'.mat'],'Prony');
        savefig([DATADIR,WaveDirection,'.fig'])
        Descriptor = fliplr(strtok(fliplr(DATADIR),'\'));
        Prony.ExportToExcel([DATADIR,WaveDirection,'.xls'],['Directory: ',Descriptor,' WaveDirection: ',WaveDirection]) 
    end
end