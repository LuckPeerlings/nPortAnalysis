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
    case 'ComsolAir'
        T = 293.15;
        R = 287;
        Properties.AmbPressure = 101.325e3;
        Properties.AmbTemperature = T-273.15;
        Properties.Viscosity = -8.38278E-7+8.35717342E-8*T^1-7.69429583E-11*T^2+4.6437266E-14*T^3-1.06585607E-17*T^4;
        Properties.SpcfHeatCap = 1047.63657-0.372589265*T^1+9.45304214E-4*T^2-6.02409443E-7*T^3+1.2858961E-10*T^4;
        Properties.Density = P0*0.02897/8.314/T;
        Properties.ThermCond = -0.00227583562+1.15480022E-4*T^1-7.90252856E-8*T^2+4.11702505E-11*T^3-7.43864331E-15*T^4;
        Properties.SpeedOfSound = sqrt(1.4*R*T);
        Properties.Xc = 0;
        Properties.Prandtl = Properties.Viscosity .* Properties.SpcfHeatCap./Properties.ThermCond;
        Properties.Gamma = Properties.SpcfHeatCap/ ( Properties.SpcfHeatCap - R);
        Properties.AmbHumidity = 0;
        Properties.ThermDiff = Properties.ThermCond./(Properties.Density .* Properties.SpcfHeatCap);
    otherwise
        error('The string in GasName is not defined in this function');
end
          
         
end