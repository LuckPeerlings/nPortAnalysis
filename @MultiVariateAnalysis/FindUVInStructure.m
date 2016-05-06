function [LocationName] = FindFieldInStructure(Input,TestName,varargin)
%FUNCTION_NAME - One line description of what the function or script performs (H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

% Author: FirstName FamilyName
% Work address
% email: 
% Website: http://www.
% May 2004; Last revision: 12-May-2004

%------------- BEGIN CODE --------------

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

% At the current positiopn in the structure (NestedName), the field names
% are tested if the wanted value is found (TestName). As there is only one
% unique field name, the sum of the strcmp will always return 0 or 1. If
% the field name is found, the NestedName is appended with the TestName and
% stored on the LocationName cell.
if sum(strcmp(fieldnames(Input),TestName))
   LocationName{length(LocationName)+1} = [NestedName, TestName];
end

%Loop over the fieldnames, and if one of the field names is a structure, go
%in to this structure and call this function again. The position is
%remembered in NestedName which is updated when the function is called. The
%base position, the name of the current structure is saved and when a new
%field name is evaluated, the NestedName is reset.
 fn = fieldnames(Input);
 for ii = 1:length(fn)    
    if ii == 1; BasePosition = NestedName; end %Save the base position
    if isstruct(Input.(fn{ii})) %If a member of the structure is a structure, save the name of that field and go deeper into the structure
    NestedName{length(NestedName)+1} = fn{ii};    
    [LocationName] = fnc_FindFieldInStructure(Input.(fn{ii}),TestName,LocationName,NestedName);
    end
    NestedName = BasePosition; %Reset the NestedName to the base position
 end
end

%------------- END OF CODE --------------
%Please send suggestions for improvement of the above template header 
%to Denis Gilbert at this email address: gilbertd@dfo-mpo.gc.ca.
%Your contribution towards improving this template will be acknowledged in
%the "Changes" section of the TEMPLATE_HEADER web page on the Matlab
%Central File Exchange


