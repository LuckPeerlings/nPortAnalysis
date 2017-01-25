function obj = SimulateVariance(obj)

%Determine the mean values of the output
BaseOutput = UpdateSolution(obj,obj.BaseInput);
for nn = 1:length(obj.UVOutputList)   
    Pos = obj.UVOutputList(nn).Pos;
    if isempty(obj.UVOutputList(nn).StructPos)
        SingleOutputValue = BaseOutput{Pos };
    else 
        SingleOutputValue = getfield(BaseOutput{Pos}, obj.UVOutputList(nn).StructPos{:});
    end
    obj.UVOutputList(nn).UV.Value = SingleOutputValue;  
end

%Perform the Cholesky decomposition of the covariance matrix and safe it.
%This to iterate faster when performing the MonteCarlo simulation
L = GetCholeskyDecomposition(obj);

for II = 1:obj.Iterations  
    PerturbedInput = CreatePerturbedValues(obj,L);
    PerturbedOutput = UpdateSolution(obj,PerturbedInput);
    obj = UpdateVariance(obj,PerturbedOutput,obj.TotalIterations+II);
end
obj.TotalIterations = obj.TotalIterations + obj.Iterations;
end


function L = GetCholeskyDecomposition(obj)
for nn = 1:length(obj.UVInputList)   
    UV = obj.UVInputList(nn).UV; %The uncertain variable
    for jj = 1:size(UV.Value,1)
        for ff = 1:size(UV.Value,2)
            %Determine the covariance matrix from the data
            CovarianceMatrix =  reshape(UV.UCMatrix(jj,ff,:),2,2)*...
                                reshape(UV.CorrelationMatrix(jj,ff,:),2,2)*...
                                reshape(UV.UCMatrix(jj,ff,:),2,2);
            %Perform the cholesky decomposition to create the correlation.
            %Exception when the variable isreal, as the cholesky
            %decomposition is not defined
            if CovarianceMatrix(2,2) == 0
                L{nn}{jj}{ff} = zeros(2,2);
                L{nn}{jj}{ff}(1,1) = sqrt(CovarianceMatrix(1,1));
            else
                L{nn}{jj}{ff} = chol(CovarianceMatrix,'lower');
            end
        end
    end
end
end

function PerturbedInput = CreatePerturbedValues(obj,L)
%Perturb all the values at once
PerturbedInput = obj.BaseInput;

for nn = 1:length(obj.UVInputList)   
    UV = obj.UVInputList(nn).UV; %The uncertain variable
    ii = obj.UVInputList(nn).Pos;

    for jj = 1:size(UV.Value,1)
        for ff = 1:size(UV.Value,2)
            %Create The Perturbed Value                          
            RandomVec = randn(2,1);
            PerturbedValue(jj,ff) = [1 1i]*L{nn}{jj}{ff}*RandomVec + UV.Value(jj,ff);
        end    
        %Update the perturbed input vector
        if isempty(obj.UVInputList(nn).StructPos)                    
            PerturbedInput{ii}{2} = PerturbedValue;
        else
            PerturbedInput{ii}{2} = setfield(PerturbedInput{ii}{2},obj.UVInputList(nn).StructPos{:},PerturbedValue);   
        end   
                  
    end
    %Clear the 
       PerturbedValue = [];  
end
end

function obj = UpdateVariance(obj, PerturbedOutput,IterationNumber)
%The variance is calculated
%For the second iteration 
for nn = 1:length(obj.UVOutputList)   
    Pos = obj.UVOutputList(nn).Pos;
    if isempty(obj.UVOutputList(nn).StructPos)
        NewData = PerturbedOutput{obj.UVOutputList(nn).Pos };
    else 
        NewData = getfield(PerturbedOutput{Pos}, obj.UVOutputList(nn).StructPos{:});
    end

    if IterationNumber == 1;   
        obj.UVOutputList(nn).UV.MeanValue = NewData;
        obj.UVOutputList(nn).UV.Var(:,:,:) = zeros([size(NewData),4]);    
        obj.UVOutputList(nn).UV.Value_Iteration(IterationNumber,:) = NewData;
        continue
    else
        PrevMeanValue = obj.UVOutputList(nn).UV.MeanValue;
        PrevRealVariance = reshape(obj.UVOutputList(nn).UV.Var(1,:,1),1,[]);
        PrevCoVariance = reshape(obj.UVOutputList(nn).UV.Var(1,:,3),1,[]);        
        PrevImagVariance = reshape(obj.UVOutputList(nn).UV.Var(1,:,4),1,[]);        
        
        MeanValue = (NewData + PrevMeanValue*(IterationNumber-1))/IterationNumber;
        obj.UVOutputList(nn).UV.MeanValue = MeanValue;
        obj.UVOutputList(nn).UV.Value_Iteration(IterationNumber,:) = NewData;
    end
    
    %An interative technique is used to update the mean value and the
    %correlation between the real and imaginary part
    %Update Mean Value
   
    
%     disp('sizes')
%     size(PrevMeanValue)
%     size(PrevCoVariance)
%     size(PrevRealVariance)
%     size(PrevImagVariance)
%     size(MeanValue)
    %Update real and imaginary variance and the covariance
    %For the first two iterations, the quantity M2n is saved instead of the
    %unbiased estimate of the variance
    if IterationNumber < 2  
        RealVariance = PrevRealVariance + real(NewData - PrevMeanValue).*real(NewData-MeanValue);
        ImagVariance = PrevImagVariance + imag(NewData - PrevMeanValue).*imag(NewData-MeanValue);
        CoVariance = PrevCoVariance + (IterationNumber-1)/IterationNumber *...
                                                real(NewData-PrevMeanValue).*imag(NewData-PrevMeanValue);   
    elseif IterationNumber == 2         
        RealVariance = 1/(IterationNumber-1)*(PrevRealVariance + real(NewData - PrevMeanValue).*real(NewData-MeanValue));
        ImagVariance = 1/(IterationNumber-1)*(PrevImagVariance + imag(NewData - PrevMeanValue).*imag(NewData-MeanValue));
        CoVariance =    1/(IterationNumber-1)*( PrevCoVariance + (IterationNumber-1)/IterationNumber *...
                                                                    real(NewData-PrevMeanValue).*imag(NewData-PrevMeanValue));
    else
        RealVariance = 1/(IterationNumber-1)*(PrevRealVariance*(IterationNumber-2) + real(NewData - PrevMeanValue).*real(NewData-MeanValue));
        ImagVariance = 1/(IterationNumber-1)*(PrevImagVariance*(IterationNumber-2) + imag(NewData - PrevMeanValue).*imag(NewData-MeanValue));
        CoVariance =   1/(IterationNumber-1)*( PrevCoVariance*(IterationNumber-2) + (IterationNumber-1)/IterationNumber *...
                                                            real(NewData-PrevMeanValue).*imag(NewData-PrevMeanValue));
    end
%     disp('sizes2')
%     size(CoVariance)
%     size(RealVariance)
%     size(ImagVariance)
%     size(MeanValue) 
    
    %Save the new variances to the uncertainVariables
    Var(:,1) = RealVariance;
    Var(:,2) = CoVariance;
    Var(:,3) = CoVariance;
    Var(:,4) = ImagVariance;
    obj.UVOutputList(nn).UV.Var(1,:,:) = Var;   
    % Calculate the rotated covariance matrix
    for mm = 1:size(obj.UVOutputList(nn).UV.Value,2)
        Theta = angle(obj.UVOutputList(nn).UV.Value(1,mm));
        R = [cos(-Theta), -sin(-Theta); sin(-Theta), cos(-Theta)];
        AlignedUCMatrix = R*reshape(Var(mm,:),2,2)*transp(R);
        obj.UVOutputList(nn).UV.AlignedVar(1,mm,1:4) = AlignedUCMatrix(:);
    end
    
end
            
end

