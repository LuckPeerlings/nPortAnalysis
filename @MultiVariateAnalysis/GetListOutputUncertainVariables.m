function obj = GetListOutputUncertainVariables(obj)

%BaseInput is the input for the class to the uncertainty analysis, having
%only the mean values. I.e. the UncertaintyValue class is removed and only
%the property 'Value' of that class is set.
BaseInput = obj.Input;
assignin('base','BaseInput',BaseInput);
%Determine the input structure without the uncertain variables. The Structure BaseInput contains the all the mean values of the uncertain variables  
for ii = 1:length(obj.UVInputList)    
    if isempty(obj.UVInputList(ii).StructPos)
       BaseInput{ obj.UVInputList(ii).Pos }{2}  = obj.UVInputList(ii).UV.Value;
    else        
       BaseInput{ obj.UVInputList(ii).Pos }{2} = setfield(BaseInput{ obj.UVInputList(ii).Pos }{2}, obj.UVInputList(ii).StructPos{:}...
                                                ,obj.UVInputList(ii).UV.Value); 
    end
end


obj.BaseInput = BaseInput;

%Determine the BaseOutput
BaseOutput = UpdateSolution(obj,BaseInput);
if isempty(BaseOutput)
   error('BaseOutput is empty, correct class implementation?') 
end
obj.BaseOutput = BaseOutput;

%Create the Output data structure which holds the output of the analysis.
%The Output data structure is created such that all numbers are saved as
%uncertain variables
%If the property is a direct value, then the cell array will have that
%value, Output{ii} = value, otherwise the cell array will contain the
%structure with the value properties 


% UVList{nn}{1} Gives the value of the UV in the n-th variable
% UVList{nn}{2} Gives the position of the n-th UV in the i-th index of the vector
% UVList{nn}{3} Gives the structure position of the n-th UV in the i-th index of the vector, if there is no structure this is empty

nn = 1;
for ii = 1:length(BaseOutput)
    if isstruct(BaseOutput{ii})        
        UVPos = FindValuesInStructure(BaseOutput{ii});
        for jj = 1:length(UVPos)  
            UV = UncertainVariable;
            UV.Value = getfield(BaseOutput{ii},UVPos{jj}{:});  
            if size(UV.Value,1) ~=1
                error('Output''s where the first dimension is unequal to 1 are unsupported at the moment')
            end
            UVList(nn).UV = UV;
            UVList(nn).Pos = ii;
            UVList(nn).StructPos = UVPos{jj};
            nn = nn+1;
        end        
    else
        UV = UncertainVariable;
        UV.Value = BaseOutput{ii};  
        if size(UV.Value,1) ~=1
            error('Output''s where the first dimension is unequal to 1 are unsupported at the moment')
        end
        UVList(nn).UV = UV;
        UVList(nn).Pos = ii;
        UVList(nn).StructPos = [];
        nn = nn+1;    
    end    
end

obj.UVOutputList = UVList;
end

function [LocationName] = FindValuesInStructure(Input, varargin)
% This function is called recursively and stores the information of the
% current position in the structure (NestedName), the found locations of
% the needed fieldname (LocationName). If the function is called the first
% time, theses values should not be given.
if isempty(varargin)
   LocationName = cell.empty;
   NestedName = cell.empty;
else
   LocationName = varargin{1};
   NestedName = varargin{2};
end
% At the current position in the structure (NestedName), the field names
% are tested if the wanted value is found (TestName). As there is only one
% unique field name, the sum of the strcmp will always return 0 or 1. If
% the field name is found, the NestedName is appended with the TestName and
% stored on the LocationName cell.
FieldNames = fields(Input);
for ii = 1:length(FieldNames)
    if ~isstruct(Input.(FieldNames{ii}))
      LocationName{length(LocationName)+1} = [NestedName, FieldNames{ii}];
    end
end
%Loop over the fieldnames, and if one of the field names is a structure, go
%in to this structure and call this function again. The position is
%remembered in NestedName which is updated when the function is called. The
%base position, the name of the current structure is saved and when a new
%field name is evaluated, the NestedName is reset.
fn = fields(Input);
for ii = 1:length(fn)    
    if ii == 1; BasePosition = NestedName; end %Save the base position
    if isstruct(Input.(fn{ii})) %If a member of the structure is a structure, save the name of that field and go deeper into the structure
    NestedName{length(NestedName)+1} = fn{ii};    
    [LocationName] = MultiVariateAnalysis.FindValuesInStructure(Input.(fn{ii}),LocationName,NestedName);
    end
    NestedName = BasePosition; %Reset the NestedName to the base position
end 
end