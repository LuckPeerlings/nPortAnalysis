function AirProperties_chk(varargin)
%Function to check te input to the AirProperties 

pars = inputParser;   %Create a parser object
DEFAULT.p = [];   %Default Pressure
DEFAULT.t = [];       %Default temperature   
DEFAULT.h = [];        %Default relative humidity
DEFAULT.x_c = []; %Default content of CO2 in air.
DEFAULT.GetOutput = false;  %Default value of the GetOutput flag

% The available arguments and their check function
addParameter(pars,'t',DEFAULT.t,@isfloat)
addParameter(pars,'p',DEFAULT.p,@isfloat)
addParameter(pars,'RH',DEFAULT.h,@isfloat)
addParameter(pars,'Xc',DEFAULT.x_c,@isfloat)
addParameter(pars,'GetOutput',DEFAULT.GetOutput,@islogical)

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

if isempty(Input.p); error('No pressure information is given'); end
if isempty(Input.t); error('No temperature information is given'); end
if isempty(Input.h); error('No humidity information is given'); end
if isempty(Input.x_c); error('No CO2 concentration information is given'); end


%Check wether if the values contain an vector
%Loop over the variables that are used in the computation to see if each of
%these variables is defined as a row vector and save the sizes if the
%inputs are a vector.


SizeVector = [];
fn = fieldnames(Input);
for jj = 1:length(fn)
    if strcmp(fn{jj},'GasName')
        continue
    end
    if ~isrow(Input.(fn{jj}));
        error('The input %s is not defined as a row-vector',fn{jj})
    end
    %Save the size of the input if it is a row vector
    if size(Input.(fn{jj}),2) ~= 1;
        SizeVector(length(SizeVector)+1) = size(Input.(fn{jj}),2);
    end
end
%Check if all the vectors have the same length
if isempty(SizeVector); SizeVector = 1; end
if length(unique(SizeVector))~= 1
   error('Input Arrays for the AirProperties do not have the same length')
end

%Loop over the variables, and the ones that are scalars are saved as
%vectors
for jj = 1:length(fn)
   if size(Input.(fn{jj}),2) == 1;
        Input.(fn{jj}) = ones(1,unique(SizeVector))*Input.(fn{jj});
    end
end

% Check if the inputs are correct
if sum( (Input.p < 0) ) ~= 0
    error('The absolute pressure can not be lower then 0 Pa');
end
if sum( (Input.h > 1) + (Input.h < 0)) ~= 0
    error('The humidity is not given as a number beteween 0 and 1');
end
if sum( (Input.t < -273.15) ) ~= 0
    error('The temperature can not be lower then 273.15 Celsius');
end
if sum( (Input.x_c < 0) + (Input.x_c > 1) ) ~= 0
    error('The CO2 Molar fraction should be between 0 and 1');
end

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

% Check if the correlations are valid
for ii = 1:unique(SizeVector)
if ~( (0 <= x_w(ii)) && (x_w(ii)<=0.06) )
    warning('Correlation is not valid for this mole fraction of water vapor in the air, 0 < x_w < 0.06')
end
if ~( (75000<=Input.p(ii)) && (Input.p(ii)<=102000) )
    warning('Correlation is not valid for this pressure, 75000 < p < 102000 Pascal')
end
if ~( (0<=Input.t(ii)) && (Input.t(ii)<=30) )
    warning('Correlation is not valid for this temperature, 0 < t < 30 Celsius')
end
end


end


