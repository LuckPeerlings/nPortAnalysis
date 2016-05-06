clear all
close all
addpath('../../');
TestObject = TestClass;


TestCalc = MultiVariateAnalysis;
A.x = UncertainVariable(rand(3),[1;1;1],[1;1;1]);
A.y = UncertainVariable(rand(3),[2;2;2],[1;1;1]);

A.x.setCorrelation(1,'A',3,[],'y')
TestCalc.Input{1}{1} = 'A';
TestCalc.Input{1}{2} = A;

TestCalc.OutputProperties{1} = 'z';
TestCalc.ClassHandle = TestObject;
TestCalc.MethodHandles = {'Calculate'};
TestCalc.GetListInputUncertainVariables;
TestCalc.GetListOutputUncertainVariables;
TestCalc.DetermineSensitivity;
TestCalc.UnCorrVariance;
TestCalc.CreateCorrelationInfo;
TestCalc.CorrVariance;



