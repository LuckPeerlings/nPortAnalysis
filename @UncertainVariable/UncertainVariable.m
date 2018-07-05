classdef UncertainVariable < handle
    properties (SetAccess = public)
        Value       %The mean value of the measurand
        MeanValue   %The mean value obtained from the Monte Carlo
        UCMatrix    %The (co) variance (matrix) for a real or oomplex matrix
        AlignedUCMatrix
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
%         Type = 'Input'
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
%         obj.Type = 'Input';
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
        function obj = calculateTotalUncertainty(obj)
            if strcmp(obj.Type,'Input')
                error('Cannot calculate the covariance matrix for uncertain variables of the input type');                
            end
            obj.AlignedUCMatrix = zeros(length(obj.Value),4);
            obj.UCMatrix = zeros(length(obj.Value),4);            
            for ii = 1:size(obj.Var,1)
                obj.UCMatrix = obj.UCMatrix + squeeze(obj.Var(ii,:,:));       
                obj.AlignedUCMatrix = obj.AlignedUCMatrix + squeeze(obj.AlignedVar(ii,:,:));   
            end
        
        end
        function obj = plotRelativeUncertainty_Input(obj)
            figure
            for ii = 1:size(obj.UCMatrix,1)
                AX{ii} = subplot(size(obj.UCMatrix,1),1,ii);
                set(AX{ii},'xtick',[])
                set(AX{ii},'xticklabel',[])
                for jj = 1:size(obj.UCMatrix,2)
                    CovarianceMatrix = reshape(squeeze(obj.UCMatrix(ii,jj,:)),2,2);
                    [~,D] = eig(CovarianceMatrix,'vector');
                    AreaEllipse(jj) = sqrt(D(1)*D(2));
                end
                %set(gca,'xtick',[])
                %set(gca,'xticklabel',[])
                    semilogy(AreaEllipse(jj)./abs(obj.Value(ii,:)))
                    grid on
            end
        end 
        function obj = plotGroupContributionRealImag(obj,XValues)
            if nargin == 1
                XValues = [1:size(obj.Value,2)].';
            end
            GroupNames = obj.Group;
                     
            nn = 1;
            while ~isempty(GroupNames)
                GroupNameUnique(nn) = GroupNames{1};
                GroupNames(1) = [];
                
                %Remove all group names that are the same
                mm = 1;
                GroupNames_ToDelete = [];
                for jj = 1:length(GroupNames)
                    if strcmp(GroupNames{jj},GroupNameUnique{nn})
                        GroupNames_ToDelete(mm) = jj;
                        mm = mm +1;
                    end
                end
                
                GroupNames(GroupNames_ToDelete) = [];
                nn = nn + 1;
                %See which variances belong to the specific groupnames
            end
           
            CoVar_Group = zeros(length(GroupNameUnique),size(obj.Value,2),4);
             for nn = 1:length(GroupNameUnique)
                for jj = 1:length(obj.Group)
                    if strcmp(obj.Group{jj},GroupNameUnique{nn})
                        Covar_Group(nn,:,:) = CoVar_Group(nn,:,:) + obj.Var(jj,:,:);                        
                    end
                    
                end
             end
             figure;
             subplot(1,2,1)
             area(XValues, squeeze(Covar_Group(:,:,1).'))
             legend( GroupNameUnique{:})                      
             title('Variance distribution Real Value')
             
             subplot(1,2,2)
             area(XValues, squeeze(Covar_Group(:,:,4).'))
             legend( GroupNameUnique{:})               
             title('Variance distribution Imaginary Value')
        end
        function [AX1,AX2] = plotRealImag(obj,XValues,SD) 
           
            Value = obj.Value.';
            if isempty(obj.CorrVar)
                TotVar = squeeze(sum(obj.Var,1));
            else
                TotVar = squeeze(sum(obj.Var,1) + sum(obj.CorrVar,1));
            end
            if nargin == 1
                XValues = [1:length(Value)].';
                SD = 1.65;
            end
            if nargin == 2
                SD = 1.65;
            end
            %Calculate the variance matrix rotated in the direction of
            %the mean vector value.
            %Taken from In-Phase/Quadrature Covariance-Matrix
            %Representation of the Uncertainty of Vectors and Complex
            %Numbers, Dylan F. Willians, C.m. Wand and Uwe Arz
            figure;
            Top = real( Value) + SD*sqrt(TotVar(:,4));
            Bottom = flipud(real( Value) - SD*sqrt(TotVar(:,4)));


            AX1 = subplot(1,2,1);
            hold all            
            fill([XValues;flipud(XValues)],[Top;Bottom],[0.8,0.8,0.8],'EdgeColor','none')
            plot(XValues,real( Value),'b-');
        
            
            Top = imag( Value) + SD*sqrt(TotVar(:,4));
            Bottom = flipud(imag( Value) - SD*sqrt(TotVar(:,4)));
            
            AX2 = subplot(1,2,2);
            hold all            
            fill([XValues;flipud(XValues)],[Top;Bottom],[0.8,0.8,0.8],'EdgeColor','none')
            plot(XValues,imag( Value),'b-');
        
        end
        
        function obj = plotRelativeUncertainty_Output(obj)
            figure
            for ii = 1:size(obj.Var,1)
                AX{ii} = subplot(size(obj.Var,1),1,ii);
                set(AX{ii},'xtick',[])
                set(AX{ii},'xticklabel',[])
                for jj = 1:size(obj.Var,2)
                    CovarianceMatrix = reshape(squeeze(obj.Var(ii,jj,:)),2,2);
                    [~,D] = eig(CovarianceMatrix,'vector');
                    AreaEllipse(jj) = sqrt(D(1)*D(2));
                end
                %set(gca,'xtick',[])
                %set(gca,'xticklabel',[])
                    semilogy(AreaEllipse(jj)./abs(obj.Value(ii,:)))
                    grid on
            end
        end 
        function obj = plotUncertaintyEllipse(obj,axesHandle,SD,Color)
            if isreal(obj.Value)
                error('UncertaintyEllipse can not be plotted for real numbers')
            end
            if nargin == 3
                Color = [0.6,0.6,0.6];
            end
            if nargin == 1
                axeshandle = gca;
                SD = 2;
                Color = [0,0,0];
            end
            for ii=1:size(obj.Value,2)
            H = fnc_plot_gaussian_ellipsoid(   [real(obj.Value(:,ii)),imag(obj.Value(:,ii))],...
                                            reshape(obj.UCMatrix(ii,:),[2,2]),SD,100,axesHandle);
                                        set(H,'color',Color);
            end
        %  plots the distribution specified by 
        %  mean M and covariance C. The distribution is plotted as an ellipse (in 
        %  2-d) or an ellipsoid (in 3-d).  By default, the distributions are 
        %  plotted in the current axes. H is the graphics handle to the plotted 
        %  ellipse or ellipsoid.
        end       
        function removeCorrelation(Position)
        end
        
        %Overloaded operators
        %The operators are overloaded such that it functions exactly like a
        %double.
        %The list is not exhaustive and operators have to be added
        %depending on the need of the current classes.
        function r = sqrt(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = builtin('sqrt',a);                     
        end
        function r = exp(obj1)
            r = UncertainVariable;
            r.Value = builtin('exp',obj1.Value);                     
        end
        function r = tan(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = builtin('tan',a);          
        end
        function r = real(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = builtin('real',a);          
        end
        function r = imag(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = builtin('imag',a);          
        end
        function plot(varargin)
                     
            if nargin == 1
                [a] = CheckOverload(varargin{1},[]);
                plot(a);
            elseif nargin == 2
                if isa(varargin{1},'Axes')
                    [a] = CheckOverload(varargin{2},[]);
                    plot(varargin{1},a);
                else
                    [a,b] = CheckOverload(varargin{1},varargin{2});
                    plot(a,b);
                end
            elseif nargin == 3
                if isa(varargin{1},'Axes')
                    [a,b] = CheckOverload(varargin{2},varargin{3});
                    plot(varargin{1},a,b);
                else
                    [a,b] = CheckOverload(varargin{1},varargin{2});
                    plot(a,b,varargin{3});
                end
            elseif nargin > 4
                if isa(varargin{1},'Axes')
                    [a,b] = CheckOverload(varargin{2},varargin{3});
                    plot(varargin{1},a,b,varargin{4:end});
                else
                    [a,b] = CheckOverload(varargin{1},varargin{2});
                    plot(a,b,varargin{3},varargin{4:end});
                end     
            end
        end
        function r = isrow(obj1)            
            r = builtin('isrow',obj1.Value);          
        end
        function r = max(obj1)
            r = UncertainVariable;
            r.Value =  builtin('max',obj1.Value);          
        end
        function r = abs(obj1)
            r = UncertainVariable;
            r.Value = builtin('abs',obj1.Value);          
        end
        function r = double(obj1)
            r = builtin('double',obj1.Value);          
        end
        function r = length(obj1)
            r = builtin('length',obj1.Value);          
        end
        function r = size(obj1,varargin)            
            r = builtin('size',obj1.Value,varargin{:});          
        end
        function r = plus(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a + b;                     
        end        
        function r = minus(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a - b;                     
        end        
        function r = uminus(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r.Value = -a;                     
        end    
        function r = uplus(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = +a;                     
        end         
        function r = times(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a .* b;                     
        end
        function r = mtimes(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a * b;                     
        end        
        function r = rdivide(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a./b;                     
        end        
        function r = ldivide(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a.\b;                     
        end        
        function r = mrdivide(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a/b;                     
        end
        function r = mldivide(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a\b;                     
        end
        function r = power(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a.^b;                     
        end
        function r = mpower(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a^b;                     
        end
        function r = lt(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a<b;                     
        end
        function r = gt(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a>b;                     
        end
        function r = le(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a<=b;                     
        end
        function r = ge(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a>=b;                     
        end
        function r = ne(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a~=b;                     
        end
        function r = eq(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a==b;                     
        end
        function r = and(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a & b;                     
        end
        function r = or(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a|b;                     
        end
        function r = not(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = ~a;                     
        end
        function r = colon(obj1,d,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = a:d:b;                     
        end
        function r = ctranspose(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = a';                     
        end        
        function r = transpose(obj1)
            [a,~] = CheckOverload(obj1,[]);
            r = UncertainVariable;
            r.Value = a.';                     
        end        
        function r = horzcat(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = [a b];                     
        end
        function r = vertcat(obj1,obj2)
            [a,b] = CheckOverload(obj1,obj2);
            r = UncertainVariable;
            r.Value = [a;b];                     
        end       
        function varargout = subsref(obj,s)

            switch s(1).type
                case '.'
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                case '()'
                    if length(s) == 1
                        % Implement obj(indices)
                        r = UncertainVariable;                      
                        r.Value = builtin('subsref',obj.Value,s);
                        [varargout{1}] = r;
                    else
                        % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',obj,s);
                    end
                case '{}'
                    [varargout{1:nargout}] = builtin('subsref',obj,s);
                    
                otherwise
                    error('Not a valid indexing expression')
            end
        end
        function r = subsasgn(A,s,varargin)

            switch s(1).type
                case '.'
                    r = builtin('subsasgn',A,s,varargin{:});   
                case '()'
                    if ~isa(A,'UncertainVariable')
                        r = UncertainVariable;
                        r.Value = A;
                    else
                        r = A;
                    end
                    for ii = 1:length(varargin)
                        if isa(varargin{ii},'UncertainVariable')
                            varargin_New{ii} = varargin{ii}.Value;
                        else
                            varargin_New{ii} = varargin{ii};
                        end
                    end
                    r.Value = builtin('subsasgn',r.Value,s,varargin_New{:});                    
                case '{}'                    
                    r = builtin('subsasgn',A,s,varargin);
                otherwise
                    error('Not a valid indexing expression')
            end
        end        
        function r = subsindex(obj,s)

            r = UncertainVariable;
            r.Value = builtin('subsindex',A.Value,s);
        end        
    end
end
     

function [a,b] = CheckOverload(obj1,obj2)
    if isa(obj1,'UncertainVariable') 
        a = obj1.Value;
    else
        a = obj1;
    end
    if isempty(obj2)
        b = [];        
    else       
        if isa(obj2,'UncertainVariable') 
            b = obj2.Value;
        else
            b = obj2;
        end
    end
end
function [UCMatrix,CorrMatrix] = GetCorrUCMatrix(CovarMatrix)
     
    UCMatrix = [sqrt(CovarMatrix(1,1)), 0 , 0, sqrt(CovarMatrix(2,2))];
    
    if UCMatrix(4) == 0
        CorrMatrix = [1, 0,0,0];
    elseif UCMatrix(1) == 0
        CorrMatrix = [0, 0,0,1];
    else
        CorrMatrix = [1, CovarMatrix(1,2)/UCMatrix(1), CovarMatrix(2,1)/UCMatrix(4), 1];
    end

end