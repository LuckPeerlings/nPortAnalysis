function [LocationName] = FindObjectInStructure(Input,ObjectName,varargin)
% Function to loop over the given structure and see whether a specific
% field is present. The output is a cell array with the locations w.r.t to
% the base structure.

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
% are tested if the wanted value is found (TestName). If it is of the
% correct objectype, it's position is saved in the LocationName
% If one of the fieldnames is a structure, that field is passed again
% through this function, together with the current lists of locations and
% the current position of the field name.
FieldNames = fieldnames(Input);
for ii = 1:length(FieldNames)
    if strcmp( class(Input.(FieldNames{ii})), ObjectName)
      LocationName{length(LocationName)+1} = [NestedName, FieldNames{ii}];
    end
    if isstruct(Input.(FieldNames{ii}))
        %Create a new nestedname variable, otherwise the old one will be
        %updated after each new found structure in the input variable
        NestedNameNew = NestedName;
        NestedNameNew{length(NestedNameNew)+1} = FieldNames{ii}; 
        [LocationName] = MultiVariateAnalysis.FindObjectInStructure(Input.(FieldNames{ii}),ObjectName,LocationName,NestedNameNew);
    end
    
end

end

%------------- END OF CODE --------------
%Please send suggestions for improvement of the above template header 
%to Denis Gilbert at this email address: gilbertd@dfo-mpo.gc.ca.
%Your contribution towards improving this template will be acknowledged in
%the "Changes" section of the TEMPLATE_HEADER web page on the Matlab
%Central File Exchange


