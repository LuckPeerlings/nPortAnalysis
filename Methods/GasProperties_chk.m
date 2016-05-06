function GasProperties_chk(Input)
%Function to check if the given input is correct for the GasProperties
GasNameList = {'Air'...
              };
if ~isfield(Input,'GasName')
    error('The field name GasName is not present in the input structure')
end
          
if sum(strcmp(GasNameList, Input.GasName))~=1
  error('The string in GasName is not defined in this function');
end

switch Input.GasName
    case 'Air'
        AirProperties_chk(Input);
end

end