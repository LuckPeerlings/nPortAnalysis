function obj = CreateCorrelationInfo(obj)

%This function checks if all the correlations are set correctly

%The correlated variances are stored in the following way:
%CorrelatedVar(nn).CorrMatrix{1} Correlation matrix
%CorrelatedVar(nn).VarName Variabel name
%CorrelatedVar(nn).Pos (This is added by this script)
%CorrelatedVar(nn).Index The rows of this variabel that are being correlated
%CorrelatedVar(nn).CorrIndex The rows to which this variabel is correlated
%CorrVar(nn).StructPos If the variabel is a structure, the position of the
%               correlated variabel

for ll = 1:length(obj.Input)
    VariabelNames{ll}= obj.Input{ll}{1};
end

%Check if the indicated correlated variable is contained in the input of
%the class. Add the position of the correlated variabel to the
%CorrVar structure
for nn = 1:length(obj.UVInputList)
    if isempty(obj.UVInputList(nn).UV.CorrVar)
        %If there is not correlated variables, skip to the next variable
        continue
    end
    
    for mm = 1:length(obj.UVInputList(nn).UV.CorrVar)
        %Check if the correlated variabel exist in the list of input
        %variabels

        CheckVec = strcmp(VariabelNames,obj.UVInputList(nn).UV.CorrVar(mm).VarName);
        
        %There should be only one unique input variable name
        if sum(CheckVec)
            %Get the indices of the position in the input variable
            ii = find(CheckVec);
            obj.UVInputList(nn).UV.CorrVar(mm).Pos = ii;
            if ~isempty(obj.UVInputList(nn).UV.CorrVar(mm).StructPos)
               %Check if the structure position exists in the input
               %structure. isfield does not work with nested field,
               %therefore a work around is used with getfield
               try
                   getfield(obj.Input{ii}{2},obj.UVInputList(nn).UV.CorrVar(mm).StructPos{:});
               catch
                   error('The variable %s does not contain the indicated structure position %s',obj.UVInputList(nn).UV.CorrVar(mm).VarName,obj.UVInputList(nn).UV.CorrVar(mm).StructPos{:});
               end
               %Check if the correlation index is smaller or equal to the size of the first dimension of the input variabel is
               InputVar = getfield(obj.Input{ii}{2},obj.UVInputList(nn).UV.CorrVar(mm).StructPos{:});
               if size(InputVar,1) > obj.UVInputList(nn).UV.CorrVar(mm).CorrIndex
                   error('The indicated CorrIndex %i is not existent in %s',obj.UVInputList(nn).UV.CorrVar(mm).CorrIndex,obj.UVInputList(nn).UV.CorrVar(mm).StructPos{:});
               end
            elseif size(obj.Input{ii}{2},1) > obj.UVInputList(nn).UV.CorrVar(mm).CorrIndex
               %Check if the correlation index is smaller or equal to the size of the first dimension of the input variabel is
               error('The indicated CorrIndex %i is not existent in %s',obj.UVInputList(nn).UV.CorrVar(mm).CorrIndex,obj.UVInputList(nn).UV.CorrVar(mm).StructPos{:});
            end
        else
            error('The variable name %s is not present in the input list',obj.UVInputList(nn).UV.CorrVar(mm).VarName);
        end        
    end
end




%Create a list for which all the information is stored to find the
%positions of the correlated matrices

VarOneOutputPos = 1;
qq = 1;
for nn = 1:length(obj.UVInputList)
    %Loop over the List of Input variables and determine which index
    %belongs to the correlated variabel.
    if isempty(obj.UVInputList(nn).UV.CorrVar)
        %If there is not correlated variables, skip to the next variabel
        continue
    end    
    for mm = 1:length(obj.UVInputList(nn).UV.CorrVar)
            %Save the list position of the UV variabel, the corresponding
            %row index, the position of the correlated UV variabel and the
            %corresponding row index
            
            CorrelationInfo(qq).VarOne.InputListPos = nn;
            CorrelationInfo(qq).VarOne.Index = obj.UVInputList(nn).UV.CorrVar(mm).Index;
            VarOneOutputPos = 0;
           for ll = 1:CorrelationInfo(qq).VarOne.InputListPos-1
                VarOneOutputPos = VarOneOutputPos + size(obj.UVInputList(ll).UV.Value,1);
            end            
            CorrelationInfo(qq).VarOne.OutputListPos = VarOneOutputPos+ CorrelationInfo(qq).VarOne.Index;
            
            CorrelationInfo(qq).VarTwo.Index = obj.UVInputList(nn).UV.CorrVar(mm).CorrIndex;  
            CorrelationInfo(qq).CorrMatrix = obj.UVInputList(nn).UV.CorrVar(mm).CorrMatrix;
            
            %Find the index of the correlated variable in the UVInputList
            ii = find([obj.UVInputList(:).Pos] == obj.UVInputList(nn).UV.CorrVar(mm).Pos);
            %As the specific index could be a structure with more input
            %variables, find the matching position. There should only be
            %one such position.  
            if length(ii) > 1
                for ll = 1:length(ii)
                    if  strcmp( strjoin([obj.UVInputList(ll).StructPos],''),...
                                strjoin(obj.UVInputList(nn).UV.CorrVar(mm).StructPos,''));
                        CorrelationInfo(qq).VarTwo.InputListPos = ll;
                    end
                end
            else
                CorrelationInfo(qq).VarTwo.InputListPos = ii;
            end
            if (CorrelationInfo(qq).VarOne.InputListPos == CorrelationInfo(qq).VarTwo.InputListPos) ...
                    && (CorrelationInfo(qq).VarOne.Index == CorrelationInfo(qq).VarTwo.Index)
                error('The variabel is self correlated')
            end
            VarTwoOutputPos = 0;
            for ll = 1:CorrelationInfo(qq).VarTwo.InputListPos-1
                VarTwoOutputPos = VarTwoOutputPos + size(obj.UVInputList(ll).UV.Value,1);
            end
            CorrelationInfo(qq).VarTwo.OutputListPos = VarTwoOutputPos + CorrelationInfo(qq).VarTwo.Index;
            qq = qq + 1;
    end 
    
end
if ~exist('CorrelationInfo','var')
    warning('No correlated Variabels found')
    return
end

obj.CorrelationInfo = CorrelationInfo;
end