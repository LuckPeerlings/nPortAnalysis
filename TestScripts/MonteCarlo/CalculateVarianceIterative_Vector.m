function [TotVar,VarReal] = CalculateVarianceIterative()

load('TestMonteCarloData_Samples_8.mat')

TestData = MVA.UVOutputList.UV.Value_Iteration;

for ii = 1:size(TestData,2)
    Var(ii,:) = reshape(cov(real(TestData(:,ii)),imag(TestData(:,ii))),1,[]);
end
figure; h = plot(Var,'+')
clear Var;


MeanVal = zeros(1,size(TestData,2));
VarVal = zeros(1,size(TestData,2),4);
for ii = 1:size(TestData,1)
     [MeanVal,VarVal] = UpdateVariance(MeanVal,VarVal,TestData(ii,:),ii);
end
set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
hold all
plot(squeeze(VarVal(1,:,:)),'o')

set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
plot(squeeze(MVA.UVOutputList.UV.Var(1,:,:)),'*')



end

function [MeanVal,VarVal] = UpdateVariance(MeanVal,VarVal,NewData,IterationNumber)
%The variance is calculated
%For the second iteration 
    if IterationNumber == 1;
        MeanVal = NewData;
        VarVal(1,:,:) = zeros([size(NewData),4]);
        return
    else
        PrevMeanValue = MeanVal;
        PrevCoVariance = VarVal(1,:,3);
        PrevRealVariance = VarVal(1,:,1);
        PrevImagVariance = VarVal(1,:,4);
        MeanVal = (NewData + PrevMeanValue*(IterationNumber-1))/IterationNumber;
    end
          
    %An interative technique is used to update the mean value and the
    %correlation between the real and imaginary part
%     Update Mean Value
    
    
    %Update real and imaginary variance and the covariance
    %For the first two iterations, the quantity M2n is saved instead of the
    %unbiased estimate of the variance
    if IterationNumber < 2  
        RealVariance = PrevRealVariance + real(NewData - PrevMeanValue).*real(NewData-MeanVal);
        ImagVariance = PrevImagVariance + imag(NewData - PrevMeanValue).*imag(NewData-MeanVal);
        CoVariance = PrevCoVariance + (IterationNumber-1)/IterationNumber *...
            real(NewData-PrevMeanValue).*imag(NewData-PrevMeanValue);   
    elseif IterationNumber == 2         
        RealVariance = 1/(IterationNumber-1)*(PrevRealVariance + real(NewData - PrevMeanValue).*real(NewData-MeanVal));
        ImagVariance = 1/(IterationNumber-1)*(PrevImagVariance + imag(NewData - PrevMeanValue).*imag(NewData-MeanVal));
        CoVariance =    1/(IterationNumber-1)*( PrevCoVariance + (IterationNumber-1)/IterationNumber *...
                        real(NewData-PrevMeanValue).*imag(NewData-PrevMeanValue));
    else
        RealVariance = 1/(IterationNumber-1)*(PrevRealVariance*(IterationNumber-2) + real(NewData - PrevMeanValue).*real(NewData-MeanVal));
        ImagVariance = 1/(IterationNumber-1)*(PrevImagVariance*(IterationNumber-2) + imag(NewData - PrevMeanValue).*imag(NewData-MeanVal));
        CoVariance =   1/(IterationNumber-1)*( PrevCoVariance*(IterationNumber-2) + (IterationNumber-1)/IterationNumber *...
                        real(NewData-PrevMeanValue).*imag(NewData-PrevMeanValue));
    end
    
    %Save the new variances to the uncertainVariables
    VarVal(1,:,1) = RealVariance;
    VarVal(1,:,2) = CoVariance;
    VarVal(1,:,3) = CoVariance;
    VarVal(1,:,4) = ImagVariance;
end
