function [DecompP,Correction] = WaveDecomposition(varargin)
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
DEFAULT.CoVar = [];
DEFAULT.CompCoVar = [];
DEFAULT.f = [];        %Default frequency
DEFAULT.x = [];
DEFAULT.GasProp.t = 20;   %Default Gas Properties
DEFAULT.GasProp.GasName = 'Air';
DEFAULT.WaveNumberProp =  [];
DEFAULT.Method = 'Standard';
DEFAULT.OptimMethod = 'None'; % The optimization method that has to be used.
DEFAULT.GetOutput = false;
DEFAULT.Correction = []; %See if the correction has been calculated, if so don't 
                   %apply the optimization scheme but use the correction directly.
DEFAULT.Corr = [];

addParameter(pars,'f',DEFAULT.f)
addParameter(pars,'P',DEFAULT.P)
addParameter(pars,'CoVar',DEFAULT.CoVar)
addParameter(pars,'CompCoVar',DEFAULT.CompCoVar)
addParameter(pars,'x',DEFAULT.x)
addParameter(pars,'Method',DEFAULT.Method)
addParameter(pars,'OptimMethod',DEFAULT.OptimMethod)
addParameter(pars,'GasProp',DEFAULT.GasProp)
addParameter(pars,'WaveNumberProp',DEFAULT.WaveNumberProp)
addParameter(pars,'GetOutput',DEFAULT.GetOutput)
addParameter(pars,'Correction',DEFAULT.Correction)
addParameter(pars,'Corr',DEFAULT.Corr)

parse(pars,varargin{:});

%The pars.Results is read only, save it to a writable variable

FunctionInput = pars.Results;

%If the correction is non empty, apply the correction and skip the
%optimization
if ~isempty(FunctionInput.Correction)
    warning('The optimization method will not be used, reverting the previous determined corrections')
    switch FunctionInput.OptimMethod    
        case 'TemperatureOptimization'
            FunctionInput.GasProp.t = FunctionInput.GasProp.t + Corr.T;
        case 'TempFlowOptim'
            FunctionInput.GasProp.t = FunctionInput.GasProp.t+ Corr.T;
            FunctionInput.WaveNumberProp.U = FunctionInput.WaveNumberProp.U + Corr.U;
        otherwise
            error('Correction for this optimization procedure not defined')
    end
    FunctionInput.OptimMethod = 'None';
end
    

switch FunctionInput.OptimMethod
    case 'None'
        [DecompP,Residual] = DecompositionMethod(FunctionInput);
        Correction = [];
    case 'TemperatureOptimization'        
        %Calculate the residual with the current starting guess, find those
        %points that have a large deviation w.r.t to the standard
        %deviation of the data and don't use that data for the optimization
        x0 = 0;
        SelecVec = [];
        res = ObjectiveFunction_Temperature(x0,FunctionInput,SelecVec,false);        
        std_res = sqrt(var(res));
        SelecVec = res<std_res;       
        %In this case, the real part of the wavenumber will be optimized,
        %so the residual is as small as possible
        fun = @(x) sum(ObjectiveFunction_Temperature(x,FunctionInput,SelecVec,false));
        [x] = fminunc(fun,x0);
        %Using an unconstrained optimizer to find the right temperature.
        FunctionInput.GasProp.t = FunctionInput.GasProp.t + x;
        if abs(x)> 5;
            warning('The temperature correction is larger than 5 degrees');
            pause;
        end
        fprintf('Temperature correction %f + %f \n',mean(FunctionInput.GasProp.t),x)
        Correction.T = x;
        [DecompP,Residual] = DecompositionMethod(FunctionInput);  
   case 'TempFlowOptim'        
        %Calculate the residual with the current stating guess, find those
        %points that have a large deviation w.r.t to the standard
        %deviation of the data
        x0 = [0,0];
        SelecVec = [];
        res = ObjectiveFunction_TemperatureFlow(x0,FunctionInput,SelecVec,false);
        std_res = sqrt(var(res));
        SelecVec = res<2*std_res;
        SelecVec(1:10) = 0;
        SelecVec(end-10:end) = 0;
        %In this case, the real part of the wavenumber will be optimized,
        %so the residual is as small as possible
        fun = @(x) sum(ObjectiveFunction_TemperatureFlow(x,FunctionInput,SelecVec,false));
        
        %Using an unconstrained optimizer to find the right temperature.
        x = fminunc(fun,x0);


        if abs(x)> 5;
            warning('The temperature correction is larger than 5 degrees');
            pause;
        end
        fprintf('------\n')
        fprintf('Amount of points %i out of %i\n',sum(SelecVec),length(SelecVec))
        fprintf('Temperature correction %f + %f \n',mean(FunctionInput.GasProp.t),x(1))
        fprintf('Velocity correction %f + %f \n',FunctionInput.WaveNumberProp.U,x(2))
        Res = sum(ObjectiveFunction_TemperatureFlow(x,FunctionInput,SelecVec,false));
        Res0 = sum(ObjectiveFunction_TemperatureFlow(x0,FunctionInput,SelecVec,false));
        fprintf('Residual,Start %f, End %f \n',Res0,Res)
        FunctionInput.GasProp.t = FunctionInput.GasProp.t+x(1);
        FunctionInput.WaveNumberProp.U  = FunctionInput.WaveNumberProp.U+x(2);
        [DecompP,Residual] = DecompositionMethod(FunctionInput);        
        Correction.T = x(1);
        Correction.U = x(2);
   case 'FlowOptim'        
        %Calculate the residual with the current stating guess, find those
        %points that have a large deviation w.r.t to the standard
        %deviation of the data
        x0 = [0,0];
        SelecVec = [];
        res = ObjectiveFunction_TemperatureFlow(x0,FunctionInput,SelecVec,false);
        std_res = sqrt(var(res));
        SelecVec = res<2*std_res;
        SelecVec(1:10) = 0;
        SelecVec(end-10:end) = 0;
        %In this case, the real part of the wavenumber will be optimized,
        %so the residual is as small as possible
        fun = @(x) sum(ObjectiveFunction_Flow(x,FunctionInput,SelecVec,false));
        
        %Using an unconstrained optimizer to find the right temperature.
        x = fminunc(fun,x0);


        if abs(x)> 5;
            warning('The temperature correction is larger than 5 degrees');
            pause;
        end
        fprintf('------\n')
        fprintf('Amount of points %i out of %i\n',sum(SelecVec),length(SelecVec))
        fprintf('Velocity correction %f + %f \n',FunctionInput.WaveNumberProp.U,x(1))
        Res = sum(ObjectiveFunction_TemperatureFlow(x,FunctionInput,SelecVec,false));
        Res0 = sum(ObjectiveFunction_TemperatureFlow(x0,FunctionInput,SelecVec,false));
        fprintf('Residual,Start %f, End %f \n',Res0,Res)
        FunctionInput.WaveNumberProp.U  = FunctionInput.WaveNumberProp.U+x(1);
        [DecompP,Residual] = DecompositionMethod(FunctionInput);        
        Correction.U = x(1);        
end
end

function res = ObjectiveFunction_Temperature(x,Data,SelecVec,display)
    %The wavenumber will optimized by adding
    Data.GasProp.t = Data.GasProp.t + x;
    [DecompP,res] = DecompositionMethod(Data);
    if display
        figure; plot(res,'x'); hold all; plot(res(SelecVec),'o')
    end
    if isempty(SelecVec)
        return
    else
        res = res(SelecVec);
    end
end

function res = ObjectiveFunction_Flow(x,Data,SelecVec,display)
    %The wavenumber will optimized by adding
    Data.WaveNumberProp.U = Data.WaveNumberProp.U + x(1) ;
    [DecompP,res] = DecompositionMethod(Data);
    if display
        figure; plot(res);
    end
    if isempty(SelecVec)
        return
    else
        res = res(SelecVec);
    end
end

function res = ObjectiveFunction_TemperatureFlow(x,Data,SelecVec,display)
    %The wavenumber will optimized by adding
    Data.GasProp.t = Data.GasProp.t + x(1);    
    Data.WaveNumberProp.U = Data.WaveNumberProp.U + x(2) ;
    [DecompP,res] = DecompositionMethod(Data);
    if display
        figure; plot(res);
    end
    if isempty(SelecVec)
        return
    else
        res = res(SelecVec);
    end
end

function [DecompP,res] = DecompositionMethod(Data)
f = Data.f;
x = Data.x;
P = Data.P;
CoVar = Data.CoVar;
CompCoVar = Data.CompCoVar;
WaveNumberProp = Data.WaveNumberProp;
GasProp = Data.GasProp;
WaveNumberProp.GasProp = Data.GasProp;
WaveNumberProp.f = Data.f;

res = [];
switch Data.Method
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
    case 'Standard'
        k = NPortAnalysis.WaveNumber(WaveNumberProp);
        %Loop over the frequency vector, to set up the linear system of
        %eqations
        for ii = 1:length(f)
            for z = 1:length(x)
                A(z,:) = [exp(-1i.*k.Downstream(ii).*x(z))  exp(1i.*k.Upstream(ii).*x(z))];
            end
            b = P(:,ii);            
            %And solve either the determined or the over determined system
            Decomposition(:,ii) = (A'*A)\(A'*b);
            res(ii) = transpose(A*Decomposition(:,ii)-b)*(conj(A*Decomposition(:,ii)-b));
        end
        DecompP.Min = Decomposition(1,:);
        DecompP.Plus = Decomposition(2,:);
    case 'WLLS'
        assignin('base','WaveNumberProp',WaveNumberProp)
        k = NPortAnalysis.WaveNumber(WaveNumberProp);
        %Loop over the frequency vector, to set up the linear system of
        %eqations
        for ii = 1:length(f)
            for z = 1:length(x)
                A(z,:) = [exp(-1i.*k.Downstream(ii).*x(z))  exp(1i.*k.Upstream(ii).*x(z))];
            end
            B = [A , zeros(size(A));
                zeros(size(A)), conj(A)];
            CoVarMatrix_Aug = [squeeze(CoVar(ii,:,:)), squeeze(CompCoVar(ii,:,:)) ;
                conj(squeeze(CompCoVar(ii,:,:))), conj(squeeze(CoVar(ii,:,:)))];
            %Normalizing the augmented covariance matrix
            CoVarMatrix_Aug = CoVarMatrix_Aug/norm(CoVarMatrix_Aug);
            b_aug = [P(:,ii);conj(P(:,ii))];
            Decomposition(:,ii) = inv(B'*inv(CoVarMatrix_Aug)*B)*B'*inv(CoVarMatrix_Aug)*b_aug;
            res(ii) = transpose(B*Decomposition(:,ii)-b_aug)*conj(B*Decomposition(:,ii)-b_aug);
            
        end
        DecompP.Min = Decomposition(1,:);
        DecompP.Plus = Decomposition(2,:);
    case 'WLLS_RelWeighting'
        assignin('base','WaveNumberProp',WaveNumberProp)
        k = NPortAnalysis.WaveNumber(WaveNumberProp);
        %Loop over the frequency vector, to set up the linear system of
        %eqations
        for ii = 1:length(f)
            for z = 1:length(x)
                A(z,:) = [exp(-1i.*k.Downstream(ii).*x(z))  exp(1i.*k.Upstream(ii).*x(z))];
            end
            B = [A , zeros(size(A));
                zeros(size(A)), conj(A)];
            CoVarMatrix_Aug = [squeeze(CoVar(ii,:,:)), squeeze(CompCoVar(ii,:,:)) ;
                conj(squeeze(CompCoVar(ii,:,:))), conj(squeeze(CoVar(ii,:,:)))];
            D = norm(P(:,ii))*[ diag(1./abs(P(:,ii))), zeros(length(P(:,ii)));
                  zeros(length(P(:,ii))), diag(1./abs(P(:,ii)))];
              
            %Normalizing the augmented covariance matrix
            CoVarMatrix_Aug = D*CoVarMatrix_Aug*transpose(D)/norm(CoVarMatrix_Aug);
            b_aug = [P(:,ii);conj(P(:,ii))];
            Decomposition(:,ii) = inv(B'*inv(CoVarMatrix_Aug)*B)*B'*inv(CoVarMatrix_Aug)*b_aug;
            res(ii) = transpose(B*Decomposition(:,ii)-b_aug)*conj(B*Decomposition(:,ii)-b_aug);
            
        end
        DecompP.Min = Decomposition(1,:);
        DecompP.Plus = Decomposition(2,:);
end
if Data.GetOutput
    assignin('base','WaveDecomposition',DecompP);
end

end





