classdef UncertainVariable < handle
    properties (SetAccess = public)
        Value       %The mean value of the measurand
        MeanValue   %The mean value obtained from the Monte Carlo
        UCMatrix    %The (co) variance (matrix) for a real or oomplex matrix
        CorrelationMatrix
        CorrVar 
        DOF         %The degrees of freedom of the variable
        Identifier  %Identifier of the variable to easily track the component errors 
        Group       = {'Base'}; %Group the uncertain variable belongs to, used to group uncertainties together
        Distribution %Name of the distribution (used in monte carlo simulations)
        Value_Iteration
    end    
    properties (SetAccess = public)
        %These properties are set internally and can be obtained but not
        %set
        CrossCorrelation    %Specifies those variables which are partially correlated to this variable
        Type
    end    
    properties
        %Sensitivity to each input variable 
        Var
        AlignedVar
        Covariance
        Sensitivity
        Sensitivity_Real
        Sensitivity_Imag
        Sensitivity_InputList
        Sensitivity_Group
        d2ydx2
        error
        Index
        IndexList
        
    
    end
    methods
        %Class Constructor
        function obj = UncertainVariable(varargin)
        %The input variables are as follows
        %Varargin(1) is the mean value of the measured value, complex or real. The first dimension corresponds to input vectors, the second dimension contains information for the different frequencies.
            %		this value is mandatory.
        %Varargin(2) is the variance matrix, the first dimension correponds to the different values for a vector input, the second dimension should equal one or equal to the second dimension if the mean value.
        %		the third dimension correponds to the variance, or for a complex value the covariance matrix using column based indexing.
        %Varargin(3) correponds to the degree of freedom for each variable, the first dimension should equal to the first dimension of the mean value, second dimension equal to one or the size of the second
        %		dimension of the mean value	
        %Varargin(4) is the identifier as a cell array of strings for each
        %            variable. The size of the first dimension of the cell array should
        %            equal the first dimension of the mean value. These identifiers should
        %            be unique
        %Varargin(5) is the group indentifier, which will be used to group the end uncertainty of the measurement.
        if isempty(varargin)
            return
        end
        if length(varargin) > 6 || length(varargin) < 3
            error('The amount of input parameters is not correct')
        end  
        obj.Value = varargin{1};
        VecSize = size(obj.Value,1);        
        FreqSize = size(obj.Value,2);
            
        if ~isreal(varargin{1}); obj.Type = 'Complex'; else obj.Type = 'Real'; end
	    
        %Checks for the variance of the variable
        if  size(varargin{2},1) ~= VecSize 
            error('The length of the first dimension of the (co)variance is not equal to one or the size of the first dimension of the mean value of the UV')
        end
        if size(varargin{2},2) ~= 1 && (size(varargin{2},2) ~= FreqSize)
            error('The length of the second dimension of the (co)variance is not equal to one or the size of the second dimension of the mean value of the UV')
        end

        if strcmp(obj.Type,'Complex') && (size(varargin{2},3) ~= 4)
            %Check if the Variance has a third dimension of length 4
            %corresponding to the covariance matrix
            error('The variance of the UV does not have a third dimension of length 4 corresponding to the covariance matrix') 
        end	
        if strcmp(obj.Type,'Real') && (size(varargin{2},3) ~= 1)
            %Check if the Variance has a third dimension of length 4
            %corresponding to the covariance matrix
            error('The variance of the UV does not have a third dimension of length 1 corresponding to the variance ') 
        end
        
        %Calculate the correlation matrix and the variance matrix and save
        %it for each frequency
        if isreal(varargin{1})
            for ii = 1:size(varargin{1},1)
                if size(varargin{2},2) == 1
                    for mm = 1:FreqSize
                       [obj.UCMatrix(ii,mm,:) obj.CorrelationMatrix(ii,mm,:)] = ...
                        GetCorrUCMatrix([varargin{2}(ii,1),0;0,0]);
                    end
                else
                    for mm = 1:FreqSize
                       [obj.UCMatrix(ii,mm,:) obj.CorrelationMatrix(ii,mm,:)] = ...
                        GetCorrUCMatrix([varargin{2}(ii,mm),0;0,0]);
                    end
                    
                end
            end
        else
            for ii = 1:size(varargin{1},1)
                if size(varargin{2},2) == 1                    
                    for mm = 1:FreqSize
                    [obj.UCMatrix(ii,mm,:) obj.CorrelationMatrix(ii,mm,:)] = ...
                        GetCorrUCMatrix(reshape(varargin{2}(ii,1,:),2,2));
                    end
                else
                    for mm = 1:FreqSize
                        [obj.UCMatrix(ii,mm,:) obj.CorrelationMatrix(ii,mm,:)] = ...
                        GetCorrUCMatrix(reshape(varargin{2}(ii,mm,:),2,2));
                    end                                  
                end
            end
        end
        %Check for the DOF
        
        if  size(varargin{3},1) ~= VecSize 
            error('The length of the first dimension of the DOF is not equal to one or the size of the first dimension of the mean value of the UV')
        end
        if size(varargin{3},2) ~= 1 && size(varargin{3},2) ~= FreqSize
            error('The length of the second dimension of the DOF is not equal to one or the size of the second dimension of the mean value of the UV')
        end
        obj.DOF = varargin{3};    
        
        %Check the identifier
        if length(varargin) > 3
            if ischar(varargin{4}) 
                %Save the identifier as a cell
                obj.Identifier = varargin(4);
            elseif iscell(varargin{4}) && (length(varargin{4}) == 1 || length(varargin{4}) == VecSize)
                obj.Identifier = varargin(4);
            elseif isempty(varargin{4})
                obj.Identifier = '';
            else
                error('The fourth input is not a string or a cell and/or the length of the cell array is not equal to one of the length of the first dimension of the UV')
            end
        end
        
        %Check the group-identifier
        if length(varargin) > 4            
            if ischar(varargin{5}) 
                %Save the group identifier as a cell
                obj.Group = varargin(5);
           else
                error('The fifth input is not a string')
            end
        end
        end
        function obj = setIdentifier(obj,Indentifier)
            if ischar(Indentifier) == size(obj.Value,1)
                %Save the identifier as a cell
                obj.Identifier = Indentifier;
            elseif length(Indentifier)  == size(obj.Value,1)
                obj.Identifier = Indentifier;
            else
                error('The fourth input is not a string or a cell and/or the length of the cell array is not equal to one of the length of the first dimension of the UV')
            end
        end
        function obj = setGroup(obj,Group)
            %Function to set the group
             if ischar(Group) 
                %Save the group identifier as a cell
                obj.Group = Group;
             elseif iscell(Group)
                 if length(Group) == 1;
                     obj.Group =cellstr(Group);
                 else
                     if length(Group) == size(obj.Value,1)
                         obj.Group = cellstr(Group);
                     else
                         error('The group size is not equal to the first dimension of the uncertain variable')
                     end
                 end
             else
                 error('The input parameter is not a string')
             end            
        end

        function obj = setCorrelation(obj,CorrMatrix,VarName,Index,CorrIndex,StructPos)
            %CorrVar(nn).CorrMatrix{1} Correlation matrix
            %CorrVar(nn).VarName Variabel name
            %CorrVar(nn).Pos (This is added by this script)
            %Corrvar(nn).Index The rows of this variabel that are being correlated
            %Corrvar(nn).CorrIndex The rows to which this variabel is correlated
            %Corrvar(nn).StructPos If the variabel is a structure, the position of the
                          
            
            %The index and corr index can be a number or a vector.
            %If one of them is a number and the other is a vector, then
            %each combination the number and vector are added
            %If both of them are vectors, they should be of equal length
            %and only the combination of each vector index is made.
            %If both are a number only one combination is made.
            %If Index is a number equal to n and CorrIndex is empty, n
            %combinations are made
            if ~isempty(CorrIndex) && (length(CorrIndex) ~= length(Index))
               error('The Index and CorrIndex do not have the same length')
            end
            if isempty(CorrIndex)
                Index = [1:Index];
                CorrIndex = Index;
            end
            %The CorrMatrix can be a equal to a single number n, for this the
            %correlation matrix will be [n 0; 0; n], implying the same
            %cross correlation for both the real and imaginary part and no
            %correlation between the real and imaginary part of variabel i
            %and j
            %It can also be given as a matrix using linear indexing using the third dimension.
            %The first dimension is reserved to describe different
            %correlation matrix for the index combinations
            %The second dimension is reserved to describe different
            %correlation matrix for the independent direction
            if length(CorrMatrix) == 1
                if CorrMatrix > 1;
                    error('Invalid value for the correlation, it should be smaller than one')
                end
                CorrMatrixInput(1,1,:) = [CorrMatrix, 0, 0 ,CorrMatrix];
            else
                CorrMatrixInput = CorrMatrix;
            end
            if (size(CorrMatrixInput,2) ~= 1) && (size(CorrMatrixInput,2) ~= size(obj.Value,2))
                error('The second dimension of the CorrMatrix is not equal to that of the second dimension of the value');
            end
            if (size(CorrMatrixInput,2) ~= 1) && (size(CorrMatrixInput,3) ~= 4)
                error('The third dimension does not have the size of 4')
            end
            %The StructPos is not empty and only one string, try to split
            %it up to get the structure notation
            if (~isempty(StructPos)) && isstr(StructPos)
                StructPos = strsplit(StructPos,'.');
            end
            nn = length(obj.CorrVar) + 1;  
            
            for ii = 1:length(Index)
                if size(CorrMatrixInput,2) == 1
                    obj.CorrVar(nn).CorrMatrix = squeeze(CorrMatrixInput);
                else
                    obj.CorrVar(nn).CorrMatrix = squeeze(CorrMatrixInput(ii,:,:));
                end
                obj.CorrVar(nn).VarName = VarName;
                obj.CorrVar(nn).Index = Index(ii);
                obj.CorrVar(nn).CorrIndex = CorrIndex(ii);
                obj.CorrVar(nn).StructPos = StructPos;
                nn = nn + 1;
            end
     
        end       
        
        
        function removeCorrelation(Position)
        end
            
    end
end

function [UCMatrix,CorrMatrix] = GetCorrUCMatrix(CovarMatrix)
     
    UCMatrix = [sqrt(CovarMatrix(1,1)), 0 , 0, sqrt(CovarMatrix(2,2))];
    
    if UCMatrix(4) == 0;
        CorrMatrix = [1, 0,0,0];
    else
        CorrMatrix = [1, CovarMatrix(1,2)/UCMatrix(1), CovarMatrix(2,1)/UCMatrix(4), 1];
    end

end