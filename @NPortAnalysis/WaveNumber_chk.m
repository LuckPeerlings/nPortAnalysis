function WaveNumber_chk(varargin)
         
pars = inputParser;
DEFAULT.f = [];                 %Default frequency
DEFAULT.U =  0;                 %Default air-speed
DEFAULT.GasProp = [];           %Default gas Properties
DEFAULT.Model.Name = [];        %Default loss model.
DEFAULT.GetOutput = false;      %Default valie of the GetOutput flag

% The available arguments and their check function
addParameter(pars,'f',DEFAULT.f,@isfloat)
addParameter(pars,'U',DEFAULT.U,@isfloat)
addParameter(pars,'GasProp',DEFAULT.GasProp,@isstruct)
addParameter(pars,'Model',DEFAULT.Model)
addParameter(pars,'GetOutput',DEFAULT.GetOutput,@islogical)

% The following inputs are not used (these are the errors when using the
% sensitivity analysis and have a prefix of err_).
% One could parse all inputs, however if a variable is misspelled it will
% not be seen by the parser and not used in the calculation.


%Parsing of the arguments
parse(pars,varargin{:});
if isempty(pars.Results.f); error('The frequency vector is not defined'); end;
if ~isrow(pars.Results.f); error('The input vector for the frequency is not a row vector'); end;
if ~isrow(pars.Results.U); error('The input vector for the velocity is not a row vector'); end;
if isempty(pars.Results.f)
    error('No frequency specified')
end
if isempty(pars.Results.GasProp)
    error('The Gas Properties have not been set');
else
    % Check if an input to GasProperties is a vector, if so the length of the
    % vector should equal the length of the frequency vector.
    Prop = GasProperties(pars.Results.GasProp); 
    if size(Prop.Density,2) > 1
        if size(pars.Results.f,2) ~= size(Prop.Density,2)
            error('The length of the frequency vector is not equal to the length of the structure containing the information of the gas properties')
        end       
    end
end


%Assigning the parsed arguments to their variables

Model = pars.Results.Model;

switch Model.Name
    case 'NoLoss'           
        
    case 'Dokumaci1995' 
        %This model is taken from "Sound Transmission in Narrow Pipes with
        %Superimposed Uniform Mean Flow and Acoustic Modelling of Automobile
        %Catalytic Converters" E. Dokumaci Journal of Sound and Vibration 1995
        %vol.182 page 799-808;
        if ~isfield(Model,'r')
                error('No field name of the extra needed parameter, r, in the Dokumaci Model')
        end
        if ~isfloat(Model.r)
            error('The needed parameter ,r, is not a float')
        end        
    case 'ViscTherm2ndOrder'
        %This model is taken from "Damping and Reflection Coefficient
        %Maasurements for an Open Pipe at Low Mach and Low Helmholtz numbers"
        %MCAM Peters, A.Hirschberg, A.J. Reijnen and A.P.J. Wijnands,
        %J.Fluid.Mech 1993 vol 256 page 499-534.
        if ~isfield(Model,'r')
                error('No field name of the extra needed parameter, r, in the ViscTherm2ndOrder Model')
        end
        if ~isfloat(Model.r)
            error('The needed parameter ,r, is not a float')
        end        
    case 'Howe1995'
        if ~isfield(Model,'r')
            error('No field name of the extra needed parameter, r, in the Dokumaci Model')
        end
        if ~isfloat(Model.r)
            error('The needed parameter ,r, is not a float')
        end        
    case 'Acoubulence'
        if ~isfield(Model,'r')
            error('No field name of the extra needed parameter, r, in the Dokumaci Model')
        end
        if ~isfloat(Model.r)
            error('The needed parameter ,r, is not a float')
        end
    case 'WideDuct'
        if ~isfield(Model,'A')
            error('No field name of the extra needed parameter, A, in the wide duct Model')
        end
        if ~isfloat(Model.A)
            error('The needed parameter, A, is not a float')
        end
        if ~isfield(Model,'Perim')
            error('No field name of the extra needed parameter, Perim, in the wide duct Model')
        end
        if ~isfloat(Model.Perim)
            error('The needed parameter, Perim, is not a float')
        end
    case 'FluidLosses'
        %This model is taken from "Attenuation of sound in wide ducts with flow
        %at elevated pressure and temperature", C.Lahiri, K.Knobloch, F.Bake,
        %L.Enghardt. Journal of Sound and Vibration 2014 vol 333 page 3440-3458
        if ~isfield(Model,'r')
                error('No field name of the extra needed parameter, r, in the Dokumaci Model')
        end
        if ~isfloat(Model.r)
                error('The needed parameter ,r, is not a float')
        end        
    otherwise
        error('The string in Name is not defined in this function');
end  
   
end
