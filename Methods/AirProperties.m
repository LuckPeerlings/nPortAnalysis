function Properties = AirProperties( varargin )
%AIRPROPERTIES Properties of humid air   
%   AIRPROPERTIES gives the thermophyscial, transport properties, specific
%   heat ratio and speed of sound in humid air 
%   
%   AIRPROPERTIES() gives the properties of dry air at STP, 20 degrees celsius
%   and 314 ppm CO_2 content.
%
%   AIRPROPERTIES(...,'t',t,'p',p) gives the properties of dry air at pressure p, temperature t in degrees celsius
%   and 314 ppm CO_2 content.
%
%   AIRPROPERTIES(...'RH',h,...) gives the properties of air with a
%   relative humidity of h. h is the fraction.
%
%   AIRPROPERTIES(...'Xc',x_c,...) gives the properties of air with
%   a carbon dioxide mole fraction of x_c.
%
%   AIRPROPERTIES(...,'GetOutput',Z,...) where Z is a boolean. The output
%   of the function is saved under the variable AirProperties in the
%   workspace.
%
%   The inputs for AIRPROPERTIES can be row-vectors. If there is only
%   one row vector provided the other properties are assumed to be equal
%   for each calculated case. If there are more than one row-vector
%   provided, the row-vectors have to have the same lengths.
%
%   Output: Properties = AIRPROPERTIES
%   Properties.SpeedOfSound     Speed of sound [m/s]
%   Properties.Gamma            Specific heat ratio C_p/C_v []
%   Properties.Prandtl          Prandtl number of the mixture []
%   Properties.ThermDiff        Thermal diffusivity [m^2/s]
%   Properties.SpcfHeatCap      Specific heat capacity [J/KgK]
%   Properties.ThermCond        Thermal conductivity [W/mK]
%   Properties.Viscosity        Viscosity [Ns/m^2]
%   Properties.Density          Density [kg/m^3]
%   Properties.AmbPressure      Ambient Pressure [Pa]
%   Properties.AmbTemperature   Ambient Temperature [deg C]
%   Properties.AmbHumidity      Relative humidity [frac]
%   Properties.Xc               Concentration of CO2 [mole-frac]
%
%   The data is taken from the following references:
%   Thermophysical and transport properties of humid air at temperature range
%   between 0 and 100 degrees celsius. P.T. Tsilingiris. Energy Conversion
%   and Management 49
%   The variation of the specific heat ratio
%   and the speed of sound in air with temperature, pressure,humidity, and CO2
%   concentration. Owen Cramer. J. Acoustic. Soc. Am.93
%   
%   Author: Luck Peerlings
%   Kungliga Tekniska Högskolan, Marcus Wallenberg Laboratory, Teknikringen
%   8, 100 44 Stockholm, Sweden
%   email: luck@kth.se
%   Website: https://www.kth.se/profile/luck/
%   July 2015; Last revision: 28-July-2015
%------------- BEGIN CODE --------------

pars = inputParser;   %Create a parser object
DEFAULT.p = 101325;   %Default Pressure
DEFAULT.t = 20;       %Default temperature   
DEFAULT.h = 0;        %Default relative humidity
DEFAULT.x_c = 393.2e-9; %Default content of CO2 in air.
DEFAULT.GetOutput = false;  %Default value of the GetOutput flag

% The available arguments and their check function
addParameter(pars,'t',DEFAULT.t)
addParameter(pars,'p',DEFAULT.p)
addParameter(pars,'RH',DEFAULT.h)
addParameter(pars,'Xc',DEFAULT.x_c)
addParameter(pars,'GetOutput',DEFAULT.GetOutput)

% The following inputs are not used (these are the errors when using the
% sensitivity analysis and have a prefix of err_) and the indentifier for the specified gas (GasName).
% One could parse all inputs, however if a variable is misspelled it will
% not be seen by the parser and not used in the calcalation.
addParameter(pars,'GasName',char.empty)
%Parsing of the arguments
parse(pars,varargin{:});

%Assigning the parsed arguments to their variables
Input.p = pars.Results.p;
Input.t = pars.Results.t;
Input.h = pars.Results.RH;
Input.x_c = pars.Results.Xc;

%The next code is vectorized, if it could not be vectorized a for loop is
%used.

%%%%% First Part 
%The dependence of the specific heat ratio and speed of sound in
%air with temperature, pressure, humidity, and CO_2 concentration.

%The following data is taken form "The variation of the specific heat ratio
%and the speed of sound in air with temperature, pressure,humidity, and CO2
%concentration. Owen Cramer. J. Acoustic. Soc. Am.93 (5)

%Coefficients for the speed of sound
C(:,1) = [ 331.5024,...
             0.603055,...
            -0.000528,...
            51.471935,...
             0.1495874,...
            -0.000782,...
            -1.82e-7,...
             3.73e-8,...
            -2.93e-10,...
           -85.20931,...
            -0.228525,...
             5.91e-5,...
            -2.835149,...
            -2.15e-13,...
            29.179762,... 
             0.000486...
            ];
%Coefficients for the specific heat ratio
C(:,2) = [  1.400822,...
           -1.75e-5,...
           -1.73e-7,...
           -0.0873629,...
           -0.0001665,...
           -3.26e-6,...
            2.047e-8,...
           -1.26e-10,...
            5.939e-14,...
           -0.1199717,...
           -0.0008693,...
            1.979e-6,...
           -0.01104,...
           -3.478e-16,...
            0.0450616,...
            1.82e-6...
            ];

T = Input.t+273.15; %Thermodynamic temperature
p_sv = exp( 1.2811805e-5.*T.^2 - 1.9509874e-2*T + 34.04926034 - 6.3536311e3./T); %Saturation vapor pressure of water in air.
f = 1.00062+3.14e-8.*Input.p+5.6e-7.*Input.t.^2; %The enhancement factor
x_w = Input.h.*f.*p_sv./Input.p; %The mole fraction of water vapor in air from the relative humidity

c0    = C(1,1) + C(2,1).*Input.t + C(3,1).*Input.t.^2 + ( C(4,1) + C(5,1).*Input.t + C(6,1).*Input.t.^2 ).*x_w...
        + ( C(7,1) + C(8,1).*Input.t + C(9,1).*Input.t.^2 ).*Input.p + ( C(10,1) + C(11,1).*Input.t + C(12,1).*Input.t.^2 ).*Input.x_c...
        + C(13,1).*x_w.^2 + C(14,1).*Input.p.^2 + C(15,1).*Input.x_c.^2 + C(16,1).*x_w.*Input.p.*Input.x_c;
gamma = C(1,2) + C(2,2).*Input.t + C(3,2).*Input.t.^2 + ( C(4,2) + C(5,2).*Input.t + C(6,2).*Input.t.^2 ).*x_w...
        + ( C(7,2) + C(8,2).*Input.t + C(9,2).*Input.t.^2 ).*Input.p + ( C(10,2) + C(11,2).*Input.t + C(12,2).*Input.t.^2 ).*Input.x_c...
        + C(13,2).*x_w.^2 + C(14,2).*Input.p.^2 + C(15,2).*Input.x_c.^2 + C(16,2).*x_w.*Input.p.*Input.x_c;

%%%%% Second part

%Thermophysical and transport properties of humid air at temperature range
%between 0 and 100 degrees celsius. P.T. Tsilingiris
%Equation (7) and (8) in this pressure use the thermodynamic temperature
%which should be the normal temperature in celsius.
%Furthermore equation (9) differs significantly from the vapor pressure at
%0C celsius. Therefore it is decided to use the vapor pressure as derived
%from the previous paper

T = Input.t+273.15;   %Thermodynamic temperature [K]
P_0 = Input.p;        %Absolute pressure [Pa]
M_a = 28.9635e-3;  %Molar mass of dry air [kg/mol]
M_v = 18e-3;       %Molar mass of water [kg/mol]. Source wikipedia ;)
R = 8.3441;      %Universal gas constant [J/molK]

%%%% Thermophysical and transport properties of dry air and water vapor
mu_a = 1e-6 .* (-9.8601e-1   + 9.080125e-2 .*T    - 1.17635575e-4.*T.^2   + 1.2349703e-7  .*T.^3  - 5.7971299e-11.*T.^4); %Viscosity of dry air
k_a  = 1e0  .* (-2.276501e-3 + 1.2598485e-4.*T    - 1.4815235e-7  .*T.^2  + 1.73550646e-10.*T.^3  - 1.066657e-13 .*T.^4    + 2.47663035e-17.*T.^5);  %Thermal conductivity of dry air
c_pa = 1e3  .* ( 0.103409e1  - 0.284887e-3 .*T    + 0.7816818e-6 .*T.^2   - 0.4970786e-9  .*T.^3  + 0.1077024e-12.*T.^4)                           ; %Specific heat capacity of dry air

mu_v = 1e-7 .* ( 8.058131868e1 + 4.000549451e-1 .*Input.t);                        % Viscosity of water vapor
k_v  = 1e-3 .* ( 1.761758242e1 + 5.558941059e-2 .*Input.t + 1.663336663e-4.*Input.t.^2); % Thermal conductivity of water vapor
c_pv = 1e3  .* ( 1.86910989    - 2.578421578e-4 .*Input.t + 1.941058941e-5.*Input.t.^2); % Specific heat capacity of water vapor

P_sv = p_sv;
x_v = x_w;

A = 0.7e-8 - 0.147184e-8.*exp(1734.29./T);
B = 0.104e-14  -0.335297e-17.*exp(3645.09./T);

%%% Density of the water-vapor mixture
z_v = 1 + A.* P_sv + B.*P_sv.^2; %Compressibility factor
rho_m = 1./z_v .* P_0 ./ (R.*T) .* M_a .* (1 - x_v.*(1 - (M_v ./ M_a)));

%%% Viscosity of the water-vapor mixture
phi_av = sqrt(2)/4 .* (1+M_a./M_v).^(-0.5) .* ( 1 + (mu_a./mu_v).^0.5 * (M_v./M_a).^0.25  ).^2;
phi_va = sqrt(2)/4 .* (1+M_v./M_a).^(-0.5) .* ( 1 + (mu_v./mu_a).^0.5 * (M_a./M_v).^0.25  ).^2;
mu_m = ( (1-x_v) .* mu_a ) ./ ( (1-x_v) + x_v.*phi_av ) + ( x_v .* mu_v ) ./ ( x_v + ( 1 - x_v ).*phi_va );

%%% Thermal conductivity of the water-vapor mixture
k_m = ( (1-x_v) .* k_a ) ./ ( (1-x_v) + x_v.*phi_av ) + ( x_v .* k_v ) ./ ( x_v + ( 1 - x_v ).*phi_va );

%%% Specific heat capacity of the mixture
c_pm = ( c_pa.*(1 - x_v).*M_a + c_pv.*x_v.*M_v) ./ ( M_a.*(1 - x_v) + M_v.*x_v );

%%% Thermal diffusivity
a_m = k_m./(rho_m .* c_pm);

%%% Prandtl number
Pr_m = mu_m .* c_pm./k_m;

%%% Packaging the output in a structure
Properties.Density = rho_m;
Properties.SpeedOfSound = c0;
Properties.Viscosity = mu_m;
Properties.SpcfHeatCap = c_pm;
Properties.Gamma = gamma;
Properties.ThermCond = k_m;
Properties.Prandtl = Pr_m;
Properties.ThermDiff = a_m;
Properties.AmbPressure = Input.p;
Properties.AmbTemperature = Input.t;
Properties.AmbHumidity = Input.h;
Properties.Xc = Input.x_c;

%%% If the output is called for it will be sent to the base workspace. If
%%% the function is called multiple times, only the output of the last
%%%  call will be placed into base workspace
if pars.Results.GetOutput
    assignin('base','AirProperties',Properties);
end
end
