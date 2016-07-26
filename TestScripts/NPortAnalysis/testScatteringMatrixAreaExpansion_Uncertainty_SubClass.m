%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to calculate the scattering matrix


%close all
clear all
addpath('../');

% Loading calibration data
Z_cal  = load('C:\Users\Luck Peerlings\Documents\KTH\Experimental\MeasurementData\2013-00-Calibration\2014-11-18_8Mics.mat');

% MicPositions are taken from the rigid plate measurements, using only
% those frequencies in the plane wave range.
x_upstream =   -1*[-0.484055826894712;-0.552540158119500;-0.587711906955935;-0.757944429725209;];
x_downstream = -1*[-0.433772413098706;-0.464467905221895;-0.494905999266211;-0.680010682145556;];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the complex pressures for the upstream excitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Loading the data
    load('C:\Users\Luck Peerlings\Documents\KTH\Experimental\MeasurementData\2014-11_AreaExpansion\NoFlow\Upstream\MeasurementInfo.mat')
    load('C:\Users\Luck Peerlings\Documents\KTH\Experimental\MeasurementData\2014-11_AreaExpansion\NoFlow\Upstream\ComplexPressures.mat','MeasData');

%Applying the calibration of the measured pressures
    Z = MeasData.Z;
    for i=1:length(Measurement.Mic.Channel)
       C_fitted(:,i) = interp1(Z_cal.f_cal,Z_cal.C(:,i),MeasData.f,'Spline');
    end
    Z = Z(:,Measurement.Mic.Channel)./ C_fitted(:,Measurement.Mic.Channel);
    
% Setting the information for the wave decomposition
    
   
    Input.Port1.Constant.GasProp.Xc = 390e-9;
    Input.Port1.Constant.GasProp.GasName = 'Air';
    Input.Port1.Constant.WaveNumberProp.Model.Name = 'FluidLosses';
    Input.Port1.Constant.Method = 'Standard'; 
    Input.Port1.Constant.x = UncertainVariable(x_upstream,(x_upstream*1e-3).^2,[0;0;0;0],{'x1','x2','x3','x4'},'X1');
    Input.Port1.Constant.WaveNumberProp.Model.r = UncertainVariable(25e-3,(25e-3/100).^2,0); 

    
    Input.Port1.Meas1.GasProp.RH = Measurement.RH;
    Input.Port1.Meas1.GasProp.p = Measurement.p;    
    
    Input.Port2.Constant = Input.Port1.Constant;
    Input.Port2.Constant.WaveNumberProp.Model.r = 45e-3; 
    Input.Port2.Constant.WaveNumberProp.Model.A = pi*(45e-3)^2;    
    Input.Port2.Constant.WaveNumberProp.Model.Perim = 2*pi*45e-3;
    Input.Port2.Constant.x = UncertainVariable(x_downstream,(x_downstream*1e-3).^2,[0;0;0;0],{'x1','x2','x3','x4'},'X2');
    
    Input.Port2.Meas1.GasProp.RH = Measurement.RH;
    Input.Port2.Meas1.GasProp.p = Measurement.p;     
    
% Performing the wave decomposition for the upstream side
    Input.Port1.Meas1.GasProp.t = MeasData.T(:,1).'*100;
    Input.Port1.Meas1.P = Z(:,1:4).';
    
% Performing the wave decomposition for the downstream side
    Input.Port2.Meas1.GasProp.t = MeasData.T(:,2).'*100;
    Input.Port2.Meas1.P = Z(:,5:8).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the complex pressures for the downstream excitation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear Measurement
clear MeasData

%Loading the data
    load('C:\Users\Luck Peerlings\Documents\KTH\Experimental\MeasurementData\2014-11_AreaExpansion\NoFlow\Downstream\MeasurementInfo.mat')
    load('C:\Users\Luck Peerlings\Documents\KTH\Experimental\MeasurementData\2014-11_AreaExpansion\NoFlow\Downstream\ComplexPressures.mat','MeasData');
%Performing the calibration of the measured pressures
    Z = MeasData.Z;
    for i=1:length(Measurement.Mic.Channel)
        C_fitted(:,i) = interp1(Z_cal.f_cal,Z_cal.C(:,i),MeasData.f,'Spline');
    end
    Z = Z(:,Measurement.Mic.Channel)./ C_fitted(:,Measurement.Mic.Channel);

% Setting the information for the second measurement case which is the same
% for both ports
    
    Input.Port1.Meas2.GasProp.RH = Measurement.RH;
    Input.Port1.Meas2.GasProp.p = Measurement.p;
    
    Input.Port2.Meas2.GasProp.RH = Measurement.RH;    
    Input.Port2.Meas2.GasProp.p = Measurement.p;
    
% Setting the information for the wave decomposition for the upstream side    
    Input.Port1.Meas2.GasProp.t = MeasData.T(:,1).'*100;
    Input.Port1.Meas2.P = Z(:,1:4).';
% Performing the wave decomposition for the downstream side
    Input.Port2.Meas2.GasProp.t = MeasData.T(:,2).'*100;
    Input.Port2.Meas2.P = Z(:,5:8).';
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determining the scattering matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% oNPort = NPortAnalysis;
% oNPort.Input = Input;
% oNPort.FreqVec = MeasData.f;
% oNPort.checkInput;
% oNPort.calculateScatteringMatrix;
% oNPort.displayScatMatrix;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set the correlations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Input.Port1.Constant.x.setCorrelation(1,'Input',4,[],'Port2.Constant.x')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
oNPortUncertaintyAnalysis = UncertaintyAnalysisNPort;

oNPortUncertaintyAnalysis.Input{1}{1} = 'Input';
oNPortUncertaintyAnalysis.Input{1}{2} = Input;

oNPortUncertaintyAnalysis.Input{2}{1} = 'FreqVec';
oNPortUncertaintyAnalysis.Input{2}{2} = MeasData.f;

oNPortUncertaintyAnalysis.CalculateUncertainty;
oNPortUncertaintyAnalysis.displayScatMatrix;



