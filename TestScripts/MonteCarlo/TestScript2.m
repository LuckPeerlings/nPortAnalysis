clear all
% close all
addpath('../../');
TestObject = TestClass;


TestCalc = MonteCarlo;
% A.x = UncertainVariable([10;10;10],[1;1;1],[1;1;1]);
% A.y = UncertainVariable([10;10;10],[2;2;2],[1;1;1]);

Variance(1,1,:) = [0.5,0.02,0.02,0.1];
A.x = UncertainVariable([10+10*1i],Variance,[1]);
A.y = UncertainVariable([10],[0.0001],[1]);


%A.x.setCorrelation(1,'A',3,[],'y')
TestCalc.Input{1}{1} = 'A';
TestCalc.Input{1}{2} = A;
TestCalc.Iterations = 10;


TestCalc.OutputProperties{1} = 'z';
TestCalc.ClassHandle = TestObject;
TestCalc.MethodHandles = {'Calculate'};
TestCalc.GetListInputUncertainVariables;
TestCalc.GetListOutputUncertainVariables;

for  ii = 1:100;
TestCalc.SimulateVariance;
plot(ii,TestCalc.UVOutputList.UV.Var(1,1,1),'x')
plot(ii,TestCalc.UVOutputList.UV.Var(1,1,4),'o')
plot(ii,TestCalc.UVOutputList.UV.Var(1,1,2),'*')
hold on
end





