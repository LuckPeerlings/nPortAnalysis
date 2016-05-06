function WaveDecomposition_chk(varargin) 
% Wave Decomposition Check

pars = inputParser;
%pars.KeepUnmatched = 1; %Any input parameter is allowed
DEFAULT.P = [];
DEFAULT.f = [];        %Default frequency
DEFAULT.x = [];
DEFAULT.GasProp = [];   %Default Gas Properties
DEFAULT.WaveNumberProp =  [];    
DEFAULT.Method = 'Standard';
DEFAULT.GetOutput = false;  


addParameter(pars,'f',DEFAULT.f,@isfloat)
addParameter(pars,'P',DEFAULT.P,@isfloat)
addParameter(pars,'x',DEFAULT.x,@isfloat)
addParameter(pars,'Method',DEFAULT.Method,@isstr)
addParameter(pars,'GasProp',DEFAULT.GasProp,@isstruct)
addParameter(pars,'WaveNumberProp',DEFAULT.WaveNumberProp,@isstruct)
addParameter(pars,'GetOutput',DEFAULT.GetOutput,@islogical)

parse(pars,varargin{:});

%Check if the necessary data is available.
if isempty(pars.Results.f); error('There is no frequency data given'); end;
if isempty(pars.Results.P); error('There is no pressure data given'); end;
if isempty(pars.Results.x); error('There is no microphone position given'); end;
if isempty(pars.Results.GasProp); error('There is no information given for the gas properties'); end;
if isempty(pars.Results.WaveNumberProp); error('There is no information given to calculate the wave numbers'); end;
if ~isrow(pars.Results.f); error('The input vector for the frequency is not a row vector'); end;
if ~iscolumn(pars.Results.x); error('The input vector for the microphone positions is not a column vector'); end;

switch pars.Results.Method
    case 'Standard'
        % Perform the necessary check
        if ~iscolumn(pars.Results.x); error('The microphone position input is not a column vector'); end
        if size(pars.Results.x,1) <= 1; error('The microphone position should have at least two rows'); end
        if size(pars.Results.x,1) ~= size(pars.Results.P,1); error('The number of rows in the microphone vector is not equal to the number of rows in the pressure vector'); end
    otherwise
        error('Method name not correct')
end

end