
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

MicPositions = X.';
MicEqPositions = X.';  
P = MeasData.Z(4:end-3, SortedFrequencyIndex );
if strcmp(WaveDirection,'Upstream')
    FlowVelocity = Measurement.U*Measurement.WaveCal.Velocity_Corr.Port1;
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
