function obj = CorrVariance(obj)

if  isempty(obj.CorrelationInfo)
    warning('No correlated Variabels found')
    return
end
%Loop over all the correlations
for nn = 1:length(obj.CorrelationInfo)
    %Get the indices where the calculated sensitivties are stored
    InputPosVarOne = obj.CorrelationInfo(nn).VarOne.InputListPos;
    InputIndexVarOne = obj.CorrelationInfo(nn).VarOne.Index;
    OutputIndexVarOne = obj.CorrelationInfo(nn).VarOne.OutputListPos;
    
    InputPosVarTwo = obj.CorrelationInfo(nn).VarTwo.InputListPos;
    InputIndexVarTwo = obj.CorrelationInfo(nn).VarTwo.Index;    
    OutputIndexVarTwo = obj.CorrelationInfo(nn).VarTwo.OutputListPos;
    
    %Loop over all the outputs
    for mm = 1:length(obj.UVOutputList)
        %If the output has non singleton dimension in the second dimension
        %loop over this direction
        
        for ll = 1:size(obj.UVOutputList(mm).UV.Value,2)
            %Obtain the SensitivityMatrix for both variables
            SensitivityMatrixVarOne = reshape(obj.UVOutputList(mm).UV.Sensitivity(OutputIndexVarOne,ll,:),2,2);
            SensitivityMatrixVarTwo = reshape(obj.UVOutputList(mm).UV.Sensitivity(OutputIndexVarTwo,ll,:),2,2);
            %Obtain the correlationmatrix for the uncertain variables
            CorrelationMatrix = reshape(obj.CorrelationInfo(nn).CorrMatrix,2,2);
            %Obtain the uncertainty matrices
            if size(obj.UVInputList(InputPosVarOne).UV.UCMatrix,2) == 1
                UCMatrixOne = reshape(obj.UVInputList(InputPosVarOne).UV.UCMatrix(InputIndexVarOne,1,:),2,2);
            else
                UCMatrixOne = reshape(obj.UVInputList(InputPosVarOne).UV.UCMatrix(InputIndexVarOne,ll,:),2,2);
            end
            if size(obj.UVInputList(InputPosVarTwo).UV.UCMatrix,2) == 1
                UCMatrixTwo = reshape(obj.UVInputList(InputPosVarTwo).UV.UCMatrix(InputIndexVarTwo,1,:),2,2);
            else
                UCMatrixTwo = reshape(obj.UVInputList(InputPosVarTwo).UV.UCMatrix(InputIndexVarTwo,ll,:),2,2);
            end
            %Calculate the correlated uncertainty 
            CorrUCMatrix = SensitivityMatrixVarOne*UCMatrixOne*CorrelationMatrix*transp(SensitivityMatrixVarTwo*UCMatrixTwo);
            %The CorrUCMatrix is not diagonal!
            obj.UVOutputList(mm).UV.CorrVar(nn,ll,1:4) = CorrUCMatrix(:);
        end
    
    end
end
    