 function [DecompP] = WaveDecomposition(varargin)
% WAVEDECOMPOSITION Calculation of the Wave Decomposition
%   WAVEDECOMPOSITION calculates the the downstream and upstream
%   pressure wave.
%
%  Schematic of the definition of the geometry
%
%          P(:,1)             P(:,2)        P(:,3)     P(:,n)
%            |                  |             |          |   
% |           | <----- s ------> |             |          |
% |____________|__________________|_____________|__________|_
% |
% |      <~~~~ DecompP.Plus       
% |--> U
% |      ~~~~> DecompP.Min   
% |__________________________________________________________
%    .
%   /I\
%    I
%  x=0  x=x(1)             x=x(2)         x=x(3)        x=x(n)
%    x --> +
% x = 0 is the reference position for the decomposition and is the position
% where the object under study is located
%
% WAVEDECOMPOSITION(...,'f',x,..) gives the frequency vector for the
% decompostion. MANDOTORY
% 
% WAVEDECOMPOSITION(...,'P',x,..) gives the pressure array for the
% decompostion. MANDOTORY The amount of columns depends on the method used, but at least two
% 'TwoMicrophone':          2 Columns
% 'Overdetermination':      More than 2 columns
% 'FullWaveDecomposition':  4 Columns
%
% WAVEDECOMPOSITION(...,'x',x,..) gives the absolute position of the
% microphone with respect to the reference position for the decomposition.
% Size should equal to number of columns of the pressure array.
%
% WAVEDECOMPOSITION(...,'Method',x,..) specifies which method should be
% used to determine the up and downstream pressure waves. 
% 'Standard':
% Standard two microphone method. If length(x) = 2, two microphone method:
% Measurement of the scattering-matrix
% of acoustical two ports. Mats Åbom. Mechanical Systems and Signal
% Processing 1991 vol 5 page 89-104
% If length(x) > 2:
% On the multiple microphone method for measuring in-duct acoustic
% properties in the presence of mean flow. Seung-Ho Jang and Jeong-Guon Ih.
% J. Acoust. Soc. Am. Vol 103. 
% 'FullWaveDecomposition':  TO BE COMPLETED
% 
% WAVEDECOMPOSITION(...,'GasProp',Z,...) uses the information in Z to
% determine the speed of sound and other thermodynamic parameters of the
% acoustic media using FNC_AIRPROPERTIES. Z is a struct which contains
% the optional arguments. EX: Z.t = 20, Z.p = 101325, Z.RH = 0.5. See the
% help file of GASPROPERTIES.
%
% WAVEDECOMPOSITION(...,'WaveNumberProp',Z,...) uses the information in Z to
% determine the wavenumbers using FNC_WAVENUMBER. Z is a struct which contains
% the optional arguments. See the help file of FNC_WAVENUMBER.
%
% WAVEDECOMPOSITION(...,'GetOutput',Z,...) where Z is a boolean. When this
% function is called in another function it will parse it's output to the
% higher lying function
%
% DecompP = FNC_WAVENUMBER is the output of the function with the member
% DecompP.Upstream
% DecompP.Downstream
% If the direction of the flow is reverse, or the direction of x is
% reverse, the resulting waves will be reversed. If both are reverse then
% decompP.Upstream is still the wave travelling upstream.
%
% WARNING: PLEASE BE VERY CAREFUL WITH THE CASE SENSITIVITY OF THE INPUT
% ARGUMENTS!!!!


pars = inputParser;
DEFAULT.P = [];
DEFAULT.f = [];        %Default frequency
DEFAULT.x = [];
DEFAULT.GasProp = [];   %Default Gas Properties
DEFAULT.WaveNumberProp =  [];    
DEFAULT.Method = 'Standard';
DEFAULT.GetOutput = false;  

addParameter(pars,'f',DEFAULT.f)
addParameter(pars,'P',DEFAULT.P)
addParameter(pars,'x',DEFAULT.x)
addParameter(pars,'Method',DEFAULT.Method)
addParameter(pars,'GasProp',DEFAULT.GasProp)
addParameter(pars,'WaveNumberProp',DEFAULT.WaveNumberProp)
addParameter(pars,'GetOutput',DEFAULT.GetOutput)

parse(pars,varargin{:});

% The following inputs are not used (these are the errors when using the
% sensitivity analysis and have a prefix of err_).
% One could parse all inputs, however if a variable is misspelled it will
% not be seen by the parser and not used in the calculation.

WaveNumberProp = pars.Results.WaveNumberProp;
GasProp = pars.Results.GasProp;
WaveNumberProp.GasProp = pars.Results.GasProp;
WaveNumberProp.f = pars.Results.f;


f = pars.Results.f;
x = pars.Results.x;
P = pars.Results.P;

%PreAllocate
%Decomposition = zeros(2,length(f));

switch pars.Results.Method
    case 'Circular'
        for ii = 1:length(f)            
            ModalMatrix = NPortAnalysis.CircularModalMatrix(GasProp,WaveNumberProp,ii);
            Decomposition{ii} = ((ModalMatrix'*ModalMatrix)\(ModalMatrix'*P(:,ii)));
        end
        DecompP.Plus = zeros(length(Decomposition{end})/2,length(f));
        DecompP.Min = zeros(length(Decomposition{end})/2,length(f));
        for ii = 1:length(f) 
          %Compared to the standard definition, the p plus wave and p minus
          %wave are interchanged, such that the plus wave travels in the
          %negative x direction ( towards the measurement object)
          DecompP.Min(1:length(Decomposition{ii})/2,ii) = Decomposition{ii}(1:end/2);
          DecompP.Plus(1:length(Decomposition{ii})/2,ii) = Decomposition{ii}(end/+1:end);
        end
    case 'Rectangular'
        %Loop over the frequencies, send the data to the
        %RectangularModelMatrix to obtain the model matrices and calculate
        %the left and right running waves.
        for ii = 1:length(f)            
            ModalMatrix = NPortAnalysis.RectangularModalMatrix(GasProp,WaveNumberProp,ii);
            try
                ConditionNR(ii) = cond(ModalMatrix);
            catch
                ConditionNR(ii) = 0;
            end
            h = 0.1;
            Decomposition{ii} = (ModalMatrix'*ModalMatrix + h^2*eye(size(ModalMatrix'*ModalMatrix)))\(ModalMatrix'*P(:,ii));
        end
        DecompP.Plus = zeros(length(Decomposition{end})/2,length(f));
        DecompP.Min = zeros(length(Decomposition{end})/2,length(f));
        for ii = 1:length(f) 
          DecompP.Min(1:length(Decomposition{ii})/2,ii) = Decomposition{ii}(1:end/2);
          DecompP.Plus(1:length(Decomposition{ii})/2,ii) = Decomposition{ii}(end/2+1:end);
        end
        assignin('base','Decomp',DecompP)
        assignin('base','ConditionNR',ConditionNR)
    case 'Standard'
        k = NPortAnalysis.WaveNumber(WaveNumberProp);
        %Loop over the frequency vector, to set up the linear system of
        %eqations
        for ii = 1:length(f)
            assignin('base','WaveNumberProp',WaveNumberProp)
            for z = 1:length(WaveNumberProp.Model.x)
                A(z,:) = [exp(-1i.*k.Downstream(ii).*WaveNumberProp.Model.x(z))  exp(1i.*k.Upstream(ii).*WaveNumberProp.Model.x(z))];
            end
            b = P(:,ii);
            %And solve either the determined or the over determined system
            Decomposition(:,ii) = (A'*A)\(A'*b);
        end
        DecompP.Min = Decomposition(1,:);
        DecompP.Plus = Decomposition(2,:);   
        assignin('base','Decomp',DecompP)
end
if pars.Results.GetOutput
   assignin('base','WaveDecomposition',DecompP);
end
end
 

    

 
