%Test file for the area expansion model
%close all
clear all
AreaExpansion = AreaExpansionModel;

AreaExpansion.Temperature = 21;
AreaExpansion.Humidity= 0;
AreaExpansion.AmbientPressure = 101325;
AreaExpansion.RadiusUpstream = 25e-3;
AreaExpansion.RadiusDownstream = 45e-3;
AreaExpansion.FreqVec = linspace(100,3000,100);
AreaExpansion.ModelName = 'Kergomard';
AreaExpansion.CalculateScatteringMatrix


AreaExpansion_u = MultiVariateAnalysis;
AreaExpansion_u .Input{1}{1} = 'Temperature';
AreaExpansion_u .Input{1}{2} = UncertainVariable(21,0.1,0);

AreaExpansion_u .Input{2}{1} = 'Humidity';
AreaExpansion_u .Input{2}{2} = 0;

AreaExpansion_u .Input{3}{1} = 'AmbientPressure';
AreaExpansion_u .Input{3}{2} = 101325;

AreaExpansion_u .Input{4}{1} = 'RadiusUpstream';
AreaExpansion_u .Input{4}{2} = UncertainVariable(25e-3,(1e-5)^2,0);

AreaExpansion_u .Input{5}{1} = 'RadiusDownstream';
AreaExpansion_u .Input{5}{2} = UncertainVariable(45e-3,(5e-5)^2,0);

AreaExpansion_u .Input{6}{1} = 'FreqVec';
AreaExpansion_u .Input{6}{2} = linspace(100,2000,100);

AreaExpansion_u .Input{7}{1} = 'ModelName';
AreaExpansion_u .Input{7}{2} = 'Kergomard';


AreaExpansion_u.OutputProperties{1} = 'ScatMatrix';
AreaExpansion_u.ClassHandle = AreaExpansionModel;
AreaExpansion_u.MethodHandles = {'CalculateScatteringMatrix'};
AreaExpansion_u.GetListInputUncertainVariables;
AreaExpansion_u.GetListOutputUncertainVariables;
AreaExpansion_u.DetermineSensitivity;
AreaExpansion_u.UnCorrVariance;
AreaExpansion_u.CreateCorrelationInfo;
AreaExpansion_u.CorrVariance;

figure
for ii = 1:4
    subplot(2,2,ii)
    hold all
    plot(AreaExpansion_u .Input{6}{2}, abs(AreaExpansion_u.UVOutputList(ii).UV.Value) + 2*sqrt(squeeze(sum(AreaExpansion_u.UVOutputList(ii).UV.AlignedVar(:,:,1),1))))
    plot(AreaExpansion_u .Input{6}{2}, abs(AreaExpansion_u.UVOutputList(ii).UV.Value) + 0*sqrt(squeeze(sum(AreaExpansion_u.UVOutputList(ii).UV.AlignedVar(:,:,1),1))))
    plot(AreaExpansion_u .Input{6}{2}, abs(AreaExpansion_u.UVOutputList(ii).UV.Value) - 2*sqrt(squeeze(sum(AreaExpansion_u.UVOutputList(ii).UV.AlignedVar(:,:,1),1))))
end

figure
for ii = 1:4
    subplot(2,2,ii)
    hold all
    plot(AreaExpansion_u .Input{6}{2}, 180/pi*angle(AreaExpansion_u.UVOutputList(ii).UV.Value) + 180/pi*2*sqrt(squeeze(sum(AreaExpansion_u.UVOutputList(ii).UV.AlignedVar(:,:,4),1)))./abs(AreaExpansion_u.UVOutputList(ii).UV.Value))
    plot(AreaExpansion_u .Input{6}{2}, 180/pi*angle(AreaExpansion_u.UVOutputList(ii).UV.Value) + 180/pi*0*sqrt(squeeze(sum(AreaExpansion_u.UVOutputList(ii).UV.AlignedVar(:,:,4),1)))./abs(AreaExpansion_u.UVOutputList(ii).UV.Value))
    plot(AreaExpansion_u .Input{6}{2}, 180/pi*angle(AreaExpansion_u.UVOutputList(ii).UV.Value) - 180/pi*2*sqrt(squeeze(sum(AreaExpansion_u.UVOutputList(ii).UV.AlignedVar(:,:,4),1)))./abs(AreaExpansion_u.UVOutputList(ii).UV.Value))
end