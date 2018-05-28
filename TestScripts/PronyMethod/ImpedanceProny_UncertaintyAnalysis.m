
addpath('..\')
addpath('..\Tools')
addpath('..\Methods')

% close all
clear all

DATADIR = '.\ProcessedData\LinerC_Sample3_NoFlow_10Pa_Inserts_WideFrequency\';


WaveDirection = 'Upstream';
% WaveDirection = 'Downstream';

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
Epsilon = 25e-3;
X = 0:0.055:9*0.055;

Mic_Position = X.';
Variance_MicPosition = 1*(5e-3)^2*ones(size(X)).';
DOF_MicPosition = zeros(size(X.'));
MicPositions = UncertainVariable(  Mic_Position,...
                            Variance_MicPosition,...
                            DOF_MicPosition,...
                            [],...
                            'x'...
                            );
MicEqPositions = X.';  

 Measured_Pressure = MeasData.Z(4:end-3, SortedFrequencyIndex );
Covariance_Pressure = MeasData.CoVar_NoCorr(4:end-3,SortedFrequencyIndex ,:);
DOF_Pressure = zeros(size(Measured_Pressure,1),1);
P = UncertainVariable(  Measured_Pressure,...
                        Covariance_Pressure,...
                        DOF_Pressure,...
                        [],...
                        'P'...
                        );
% P = Measured_Pressure;
% P = Measured_Pressure;
if strcmp(WaveDirection,'Upstream')
    FlowVelocity = UncertainVariable(Measurement.U*Measurement.WaveCal.Velocity_Corr.Port1,(2)^2, 0,[],'FlowVelocity');
elseif  strcmp(WaveDirection,'Downstream')
    FlowVelocity = UncertainVariable(Measurement.U*Measurement.WaveCal.Velocity_Corr.Port2,(2)^2, 0,[],'FlowVelocity');
else
    error('WaveDirection incorrectly defined')
end

FlowVelocity = UncertainVariable(-Measurement.U*Measurement.WaveCal.Velocity_Corr.Port1,(Measurement.U*0.10)^2, 0,[],'FlowVelocity');
Temperature = UncertainVariable(Measurement.t,(1)^2, 0,[],'Temperature');
Height = UncertainVariable(0.025,(1e-3)^2, 0,[],'DuctHeight');
Prony = PronyImpedance( Measurement.f.',P,MicPositions,MicEqPositions,5,Epsilon,...
                        FlowVelocity, Temperature, Height);
                    


Prony.CalculateImpedance;
Prony.PlotImpedance;

% Z = MultiVariateAnalysis;
% Z.ClassHandle = PronyImpedance( Measurement.f.',P,MicPositions,MicEqPositions,5,Epsilon,...
%                         FlowVelocity, Temperature, Height);
% Z.MethodHandles{1} = 'CalculateImpedance';
% Z.CreateInput();
% Z.CalculateUncertainty();
% Q = Z.Output.Impedance;
% 
% [AX1,AX2] = Q.plotRealImag(Measurement.f,2);
% 
% xlabel(AX2,'frequency')

% figure; area(Q.Var(:,:,1).')
% Francois = load('C:\Users\Luck Peerlings\Documents\KTH\Experimental\MeasurementData\2018-03_LinerMeasurements_Pierre\Main folder\Impedance\Results\14-Mar-2018_LinerB_Sample1_M008_10Pa.mat');
% figure; plot(real(Francois.Input.kxProny),imag(Francois.Input.kxProny))