%Small script to merge two structures 
function Output = MergeStructs(Input1,Input2)

%The structure Input2 is going to be merged with Input 1

%Get the field names of both inputs
FN1 = fieldnames(Input1);
FN2 = fieldnames(Input2);

%Duplicate the first input to have the main structure
Output = Input1;

%Loop over the field names of input 2 and if a field name in Input1 does
%not exist it is copied, otherwise it is checked if both fields contain
%structure, then those structures are merged. Otherwise the field in Input1
%is overwritten by Input2


for ii = 1:length(FN2)
    for jj = 1:length(FN1)
        if strcmp(FN2{ii},FN1{jj}) 
            %The two structures have the same field name.
            %If both inputs are struct pass it through this same function
            %to merge those structures
            if (isstruct(Output.(FN1{jj})) && isstruct(Input2.(FN2{ii})))
                Output.(FN2{ii}) = MergeStructs(Output.(FN1{jj}),Input2.(FN2{ii}));
                Input2.(FN2{ii}) = []; %Clear the to be merged structure so that it will not be used in the next iteration
            else
                %If the fields are of different types, overwrite the main
                %structure with the information of the new field 
                if ~isempty(Input2.(FN2{ii}))
                    warning('The main structure and the to be merged structure have duplicate fields of which at least one are no structures, overwriting the field in the main structure')
                    Output.(FN2{ii}) = Input2.(FN2{ii}); 
                    Input2.(FN2{ii}) = []; %Clear the to be merged structure so that it will not be used in the next iteration
                end
            end
        else
            %If the field name does not exist, add it to the output and
            %copy the information from the other field.
            if ~isempty(Input2.(FN2{ii}))
                Output.(FN2{ii}) = Input2.(FN2{ii});
                Input2.(FN2{ii}) = []; %Clear the to be merged structure so that it will not be used in the next iteration
            end
        end
    end
end

