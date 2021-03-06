addpath('..\..\')
addpath('..\..\Methods')

close all
clear all
%% Performing the wave decomposition
%Properties of the measurement;
DATADIR = '.\Data\RigidPlate';

load([DATADIR,'\MeasurementInfo.mat'])
Source1 = load([DATADIR,'\ComplexPressures_Source1.mat'],'MeasData');


[~,Index_Meas1] = sort(Source1.MeasData.Ref_Freq);


%% Performing the calibration
Z_cal  = load([DATADIR,'\2016-08-19_RigidPlate.mat']);
for ii=1:size(Z_cal.C,2)
    Source1.C_fitted(ii,:) = polyval(Z_cal.PolReal(ii,:),Source1.MeasData.f,[],Z_cal.MU_Real(ii,:)) + ...
                            1i.*polyval(Z_cal.PolImag(ii,:),Source1.MeasData.f,[],Z_cal.MU_Imag(ii,:));
    %Remove the influence of the powerline
    Z_cal.C(Z_cal.f_cal == 50,:) = [];
    Z_cal.C(Z_cal.f_cal == 150,:) = []; 
    Z_cal.f_cal(Z_cal.f_cal == 50) = [];
    Z_cal.f_cal(Z_cal.f_cal == 150) = []; 
    
    Source1.C_fitted(ii,:) = interp1(Z_cal.f_cal, real(Z_cal.C(:,ii)),Source1.MeasData.Ref_Freq,'linear','extrap') + 1i*interp1(Z_cal.f_cal, imag(Z_cal.C(:,ii)),Source1.MeasData.Ref_Freq,'linear','extrap');
                         
end
Source1.MeasData.Z(:,:) = Source1.MeasData.Z(:,:)./ Source1.C_fitted(Measurement.Mic.Channel,:);

%% MeasurementInfo
MicPort1= [2,3,4,5,6,7,8];

Input.Port1.Constant.GasProp.Xc = 390e-9;
Input.Port1.Constant.GasProp.GasName = 'Air';
Input.Port1.Constant.WaveNumberProp.Model.Name = 'FluidLosses';
Input.Port1.Constant.WaveNumberProp.Model.r = 25e-3;
Input.Port1.Constant.Method = 'Standard'; 
Input.Port1.Constant.x = Measurement.Mic.Pos.Port1(MicPort1).';
Input.Port1.Constant.WaveNumberProp.U = 0;

Input.Port1.Meas1.GasProp.RH = Measurement.RH;
Input.Port1.Meas1.GasProp.p = Measurement.p;
Input.Port1.Meas1.GasProp.t = Measurement.t;

Input.Port1.Meas1.P = Source1.MeasData.Z(MicPort1,Index_Meas1);
             
oNPort = NPortAnalysis;
oNPort.Input = Input;
oNPort.FreqVec = Measurement.f.';
oNPort.NrPorts = 1;
oNPort.NrMeas = 1;
oNPort.calculateScatteringMatrix;
oNPort.displayScatMatrix;
