clear all
% close all
addpath('../../');
TestObject = TestClass;

%% Monte Carlo
TestCalcMC = MonteCarlo;

Variance(1,1,:) = [0.5,0.02,0.02,0.1];
A = UncertainVariable(1+1i,Variance,0);


%A.x.setCorrelation(1,'A',3,[],'y')
TestCalcMC.Input{1}{1} = 'A';
TestCalcMC.Input{1}{2} = A;
TestCalcMC.Iterations = 10;


TestCalcMC.OutputProperties{1} = 'z';
TestCalcMC.ClassHandle = TestObject;
TestCalcMC.MethodHandles = {'Calculate'};
TestCalcMC.GetListInputUncertainVariables;
TestCalcMC.GetListOutputUncertainVariables;


TestCalcMC.Iterations = 10;

% for ii = 1:100
% disp(ii)
% TestCalcMC.SimulateVariance;
% end

%% Linear variate
TestCalcMV = MultiVariateAnalysis;

TestCalcMV.Input{1}{1} = 'A';
TestCalcMV.Input{1}{2} = A;
TestCalcMV.OutputProperties{1} = 'z';
TestCalcMV.ClassHandle = TestObject;
TestCalcMV.MethodHandles = {'Calculate'};
TestCalcMV.GetListInputUncertainVariables;
TestCalcMV.GetListOutputUncertainVariables;
TestCalcMV.DetermineSensitivity;
TestCalcMV.UnCorrVariance;
TestCalcMV.SetOutput;
