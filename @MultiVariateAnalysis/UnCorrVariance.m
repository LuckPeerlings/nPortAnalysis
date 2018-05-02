function obj = UnCorrVariance(obj)


%Loop over all the output variables
for ii = 1:length(obj.UVOutputList)
    %Save the uncertain output variable for better readability
    UVOutput = obj.UVOutputList(ii).UV;
    %Loop over all the input variables
    
    ll = 1; %Reset the counter for the output variances
    for jj = 1:length(obj.UVInputList)
        %Save the uncertain input variable for better readability
        UVInput = obj.UVInputList(jj).UV;

        %Loop over the independent direction of the input variable
        for kk = 1:size(UVInput.Value,1)
            %Loop over the columns of the output (each column represent
            %an independent variable)
            for mm = 1:size(UVOutput.Sensitivity,2)
                SensitivityMatrix = reshape(UVOutput.Sensitivity(ll,mm,:),2,2);
                %Include the correlation matrix
                if size(UVInput.UCMatrix,2) == 1                    
                    UCMatrix = reshape(UVInput.UCMatrix(kk,1,:),2,2);
                    CorrelationMatrix = reshape(UVInput.CorrelationMatrix(kk,1,:),2,2);
                else    
                    UCMatrix = reshape(UVInput.UCMatrix(kk,mm,:),2,2);
                    CorrelationMatrix = reshape(UVInput.CorrelationMatrix(kk,mm,:),2,2);
                end
                %Calculate and save the covariance matrix for this specific
                %contribution
                OutputUCMatrix = SensitivityMatrix *  UCMatrix * CorrelationMatrix * transpose(SensitivityMatrix * UCMatrix);
                UVOutput.Var(ll,mm,1:4) = OutputUCMatrix(:);
                if length(UVInput.Group) == 1
                    UVOutput.Group{ll} = UVInput.Group;
                else
                    UVOutput.Group{ll} = UVInput.Group{kk};
                end
                %Calculate the variance matrix rotated in the direction of
                %the mean vector value.
                %Taken from In-Phase/Quadrature Covariance-Matrix
                %Representation of the Uncertainty of Vectors and Complex
                %Numbers, Dylan F. Willians, C.m. Wand and Uwe Arz
                Theta = angle(UVOutput.Value(1,mm));
                R = [cos(-Theta), -sin(-Theta); sin(-Theta), cos(-Theta)];
                AlignedUCMatrix = R*OutputUCMatrix*transpose(R);
                UVOutput.AlignedVar(ll,mm,1:4) = AlignedUCMatrix(:);
            end
            ll = ll + 1;
        end                
    end
    obj.UVOutputList(ii).UV = UVOutput;
end
end

