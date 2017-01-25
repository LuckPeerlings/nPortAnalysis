function obj = DetermineSensitivity(obj)
%IMPLEMENT: Add to the CalcSensitivity an extra input which gives the
%list of output values. At the moment it is calculated inside the function
%which is redundant.
%Include the correlation matrix
%IMPLEMENT: That the function can handle non structure outputs.
%IMPLEMENT: The indentifiers for the various contribiutions to the
%variance
%IMPLEMENT: FIX WHEN A VALUE IS ZERO; IT DOES NOT SUPPORT VECTORIZED
%PROBLEMS

%Values for the perturbations applied
RELPERTURBATION = 1e-6; %If the input is unequal to zero
PERTURBATION = 1e-6;    %if the input is equal to zero


%Routine to calculate the sensitivity of the output variables to each input
%variable
Index = 1;

for nn = 1:length(obj.UVInputList)     
    UV = obj.UVInputList(nn).UV; %The uncertain variable
    ii = obj.UVInputList(nn).Pos; %The index of the inputvector where the uncertain variabel is located
    for jj = 1:size(UV.Value,1)     
        PerturbedInput = obj.BaseInput;
        %Calculating the sensitivity of the output parameters do
        %perturbations in the real part
        PerturbedValue_Real = UV.Value;
        PerturbedValue_Real(jj,:) = real(UV.Value(jj,:))*(1+RELPERTURBATION) + ...
                  1i*imag((UV.Value(jj,:)));
        Perturbation = real(UV.Value(jj,:))*(RELPERTURBATION);
        if PerturbedValue_Real(jj) == 0;                    
            PerturbedValue_Real(jj) = PERTURBATION;
            Perturbation = PERTURBATION;
            disp('test')
        end
        if isempty(obj.UVInputList(nn).StructPos)                    
            PerturbedInput{ii}{2} = PerturbedValue_Real;
        else
            PerturbedInput{ii}{2} = setfield(obj.BaseInput{ii}{2},obj.UVInputList(nn).StructPos{:},PerturbedValue_Real);   
        end
        CalcSensitivity(obj,PerturbedInput,Perturbation,Index);
        %If the input is complex, calculate the sensitivity of the
        %output parameters due to perturbations in the imaginary part.
        PerturbedInput = obj.BaseInput;
        if ~isreal(UV.Value) 
            PerturbedValue_Imag = UV.Value;
            PerturbedValue_Imag(jj,:) = real(UV.Value(jj,:)) + ...
                                  1i*imag((UV.Value(jj,:)))*(1+RELPERTURBATION);
            Perturbation = RELPERTURBATION*1i*imag(UV.Value(jj,:));
            if PerturbedValue_Imag(jj) == 0;
                PerturbedValue_Imag(jj) = 1i*PERTURBATION;
                Perturbation = 1i*PERTURBATION;
            end
            if isempty(obj.UVInputList(nn).StructPos)  
                PerturbedInput = obj.BaseInput; %Reset the perturbation added by the real part
                PerturbedInput{ii}{2} = PerturbedValue_Imag;
            else
                PerturbedInput{ii}{2} = setfield(obj.BaseInput{ii}{2},obj.UVInputList(nn).StructPos{:},PerturbedValue_Imag);   
            end
            CalcSensitivity(obj,PerturbedInput,Perturbation,Index);
        end
        Index = Index + 1;
    end
end
end


function obj = CalcSensitivity(obj,PerturbedInput,Perturbation,Index)
%Calculate the perturbed output
PerturbedOutput = UpdateSolution(obj,PerturbedInput);

%Loop over the UVOutputList, which is a cell structure containing all the
%output variables.
for nn = 1:length(obj.UVOutputList)
    ii = obj.UVOutputList(nn).Pos;
    %Save the uncertain variables in each output property in a list
    if isempty(obj.UVOutputList(nn).StructPos)            
        BOSingle = obj.BaseOutput{ obj.UVOutputList(nn).Pos };
        POSingle = PerturbedOutput{ obj.UVOutputList(nn).Pos };
    else
        BOSingle = getfield(obj.BaseOutput{ii},  obj.UVOutputList(nn).StructPos{:});
        POSingle = getfield(PerturbedOutput{ii},  obj.UVOutputList(nn).StructPos{:});
    end
    %Determine if the perturbation was real or complex and calculate
    %the appropriate elements of the sensitivity matrix. Saved in
    %columnwise vector
    Sensitivity = (POSingle-BOSingle)./Perturbation;
    if isreal(Perturbation)
        obj.UVOutputList(nn).UV.Sensitivity(Index,:,1) = real(Sensitivity);
        obj.UVOutputList(nn).UV.Sensitivity(Index,:,2) = imag(Sensitivity);
        obj.UVOutputList(nn).UV.Sensitivity(Index,:,3) = 0;
        obj.UVOutputList(nn).UV.Sensitivity(Index,:,4) = 0;
    else
        obj.UVOutputList(nn).UV.Sensitivity(Index,:,3) = real(Sensitivity);
        obj.UVOutputList(nn).UV.Sensitivity(Index,:,4) = imag(Sensitivity);
    end
end
end







