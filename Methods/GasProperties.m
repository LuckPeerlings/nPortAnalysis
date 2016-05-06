function Properties = GasProperties( Input )
%GASPROPERTIES Function to obtain thermodynamic and properties for a given
%gas. 
%
% Syntax:  [Properties] = GasProperties(GasName,Input)
%
% Inputs:
%    GasName - The name of the to be used gas
%    Input - The input [structure] needed for the specific function to obtain the gas
%    properties. The field 'GasName' is mandatory and used to redirect the
%    input to a specified subfunction
%
%    Defined GasName: 'Air', redirects to AirProperties
%
% Outputs:
%    Properties - A structure containing all the properties of the gas at
%    specified conditions
%
% Example: 
%    Input.GasName = 'Air';
%    Input.t = 20;
%    [Properties] = GasProperties(Input)
%
% Other m-files required: none
% Subfunctions: AirProperties
% MAT-files required: none
%

% Author: Luck Peerlings
% Kungliga Tekniska Högskolan, Marcus Wallenberg Laboratory, Teknikringen
% 8, 100 44 Stockholm, Sweden
% email: luck@kth.se
% Website: https://www.kth.se/profile/luck/
% July 2015; Last revision: 28-July-2015

%------------- BEGIN CODE --------------
switch Input.GasName
    case 'Air'
        Properties = AirProperties(Input);
    otherwise
        error('The string in GasName is not defined in this function');
end
          
         
end