function [ WaveNumber ] = WaveNumber(varargin)
% WAVENUMBER Calculation of the wave number
%   FNC_WAVENUMBER gives the (complex) wave number
%   for various loss models for waves travelling in the downstream (plus) and
%   upstream (min) direction of a flow
%
%   WAVENUMBER(...,'f',x,..) gives the up and downstream wave number for a frequency of x [Hz] and mean
%   velocity of 0 without losses. This argument is MANDATORY and should be a scalar or row vector.
%   
%   WAVENUMBER(...,'U',x,...) gives the up and downstream wave number for a mean
%   velocity of x [m/s]. If this value is negative, the wave numbers (usptream and downstream) are
%   reversed.
%
%   WAVENUMBER(...,'Prop',Z,...) uses the information in Z to
%   determine the speed of sound and other thermodynamic parameters of the
%   acoustic media using GASPROPERTIES. Z is a struct which contains
%   the optional arguments. EX: Z.t = 20, Z.p = 101325, Z.RH = 0.5. See the
%   help file of GASPROPERTIES.
%
%   WAVENUMBER(...,'GetOutput',Z,...) where Z is a boolean. The output
%   of the function is saved under the variable WaveNumber in the
%   workspace.
%
%   WAVENUMBER(...,'Model',Q) determines the up and downstream wave
%   number with a different model. Q is a struct with the name of the model
%   and the extra needed parameters. The only model that is supported
%   includes the effect of viscous and thermal dissipation as given by the
%   paper "Sound transmission in narrow pipes with superimposed uniform
%   mean flow and acoustic modelling of automobile catalytic converters ".
%   E. Dokumaci. Journal of Sound and Vibration 1995 Vol 182 Pag 799-808.
%   To use this model use the following syntax for Q. Q.Name =
%   'Dokumaci1995', Q.r = duct radius in meter
%
%   WaveNumber  = WAVENUMBER(...) is the output of the
%   function with members WaveNumber.Upstream, WaveNumber.Downstream
%   Each column corresponds to a frequency, as of yet no higher order modes
%   are taken into account, but this can easily be achieved by expanding
%   the number of rows using linear indexing for the higher order modes.

pars = inputParser;
DEFAULT.f = [];                 %Default frequency
DEFAULT.U =  0;                 %Default air-speed
DEFAULT.PitotPressure = [];
DEFAULT.PitotSpeedCorrection = [];
DEFAULT.GasProp = [];           %Default gas Properties
DEFAULT.Model.Name = [];  %Default loss model.
DEFAULT.GetOutput = false;      %Default valie of the GetOutput flag

% The available arguments and their check function
addParameter(pars,'f',DEFAULT.f)
addParameter(pars,'U',DEFAULT.U)
addParameter(pars,'PitotPressure',DEFAULT.PitotPressure)
addParameter(pars,'PitotSpeedCorrection',DEFAULT.PitotSpeedCorrection)
addParameter(pars,'GasProp',DEFAULT.GasProp)
addParameter(pars,'Model',DEFAULT.Model)
addParameter(pars,'GetOutput',DEFAULT.GetOutput)


%Parsing of the arguments
parse(pars,varargin{:});

%Assigning the parsed arguments to their variables
Model = pars.Results.Model;
Prop = AirProperties(pars.Results.GasProp);
PitotPressure = pars.Results.PitotPressure;
PitotSpeedCorrection = pars.Results.PitotSpeedCorrection;
U = pars.Results.U;
f = pars.Results.f;
omega = 2.*pi.*f;

if ~isempty(PitotPressure) && U ~= 0
    error('Conflicting data, pressure and speed both given')
elseif ~isempty(PitotPressure)
    U = PitotPressureToVelocity(PitotPressure,PitotSpeedCorrection,Prop);
end

%Herein are all the models for the wavenumbers, each model can be accessed
%by setting the right model name, the extra needed parameters are stored in
%the fields of model
switch Model.Name
    case 'NoLoss'           
        k = omega./Prop.SpeedOfSound;
        WaveNumber.Upstream =    k./( 1-U./Prop.SpeedOfSound );
        WaveNumber.Downstream =  k./( 1+U./Prop.SpeedOfSound );
    case 'Dokumaci1995' 
        %This model is taken from "Sound Transmission in Narrow Pipes with
        %Superimposed Uniform Mean Flow and Acoustic Modelling of Automobile
        %Catalytic Converters" E. Dokumaci Journal of Sound and Vibration 1995
        %vol.182 page 799-808;
        r = Model.r;
        k = omega./Prop.SpeedOfSound;
        S = r.*sqrt(Prop.Density.*omega./Prop.Viscosity);
        K_0 = 1 + ( ( 1-1i )./( S.*sqrt(2) ) ).*( 1 + ( Prop.Gamma-1 )./( sqrt(Prop.Prandtl) ) );
        WaveNumber.Upstream = k .* K_0 ./ (1-K_0.*U./Prop.SpeedOfSound);
        WaveNumber.Downstream = k .* K_0 ./ (1+K_0.*U./Prop.SpeedOfSound);
    case 'Weng2016'
        r = Model.r;
        k = omega./Prop.SpeedOfSound;
        S = r.*sqrt(Prop.Density.*omega./Prop.Viscosity);
        alpha_0 = 1./(S.*sqrt(2)).*(1+(Prop.Gamma-1)./(sqrt(Prop.Prandtl)));
        M = U./Prop.SpeedOfSound;
        
        WaveNumber.Upstream = k.*((1+alpha_0)./(1-M) - 1i*alpha_0./(1-M).^(3/2));
        WaveNumber.Downstream = k.*((1+alpha_0)./(1+M) - 1i*alpha_0./(1+M).^(3/2));
    case 'ViscTherm2ndOrder'
        %This model is taken from "Damping and Reflection Coefficient
        %Maasurements for an Open Pipe at Low Mach and Low Helmholtz numbers"
        %MCAM Peters, A.Hirschberg, A.J. Reijnen and A.P.J. Wijnands,
        %J.Fluid.Mech 1993 vol 256 page 499-534.
        r = Model.r;
        Sh = r.*(omega .* Prop.Density./Prop.Viscosity).^(0.5);
        k = (omega./Prop.SpeedOfSound) .* ...
            (1 + (1-1i)./sqrt(2).*1./Sh.*(1 + (Prop.Gamma-1)./sqrt(Prop.Prandtl)) ...
             - 1i./Sh.^2 .* (1+(Prop.Gamma-1)./sqrt(Prop.Prandtl)-0.5.*Prop.Gamma.*(Prop.Gamma-1)./sqrt(Prop.Prandtl)));        
        WaveNumber.Upstream = k ;
        WaveNumber.Downstream = k;
    case 'Howe1995'
        freq = f;
        radius = Model.r;
        u_tau = 0.01;
        for ii = 1:length(freq)
            M0_avg  = U./Prop.SpeedOfSound(ii);
            gas.c0 = Prop.SpeedOfSound(ii);
            gas.nu = Prop.Viscosity(ii)/Prop.Density(ii);
            gas.Cp = Prop.SpcfHeatCap(ii);
            gas.gamma = Prop.Gamma(ii);
            gas.T0 = 20;
            gas.Pr = Prop.Prandtl(ii);
            gas.Pr_turb = 0.7;
            [kDown(ii), kUp(ii)] = getKFromHowe1995(freq(ii), gas, radius, M0_avg, u_tau);
        end
        WaveNumber.Upstream = kUp ;
        WaveNumber.Downstream = kDown;    
    case 'Acoubulence'         
        freqVec = f;
        radius = Model.r;
        for ff = 1:length(f)
            M0 = U./Prop.SpeedOfSound(ff);
            [kDown(ff), kUp(ff)] = fnc_Acoubulence(freqVec(ff),radius,M0);
        end
        WaveNumber.Upstream = kUp ;
        WaveNumber.Downstream = kDown;    
        
    case 'WideDuct'
        %The followin dispersion relation is obtained from Acoustics, and
        %introduction to it's physical principles and applications buy
        %Allan D. Pierce, page 534-535
        alpha_cl =  (omega).^2./ (2.*Prop.Density.*Prop.SpeedOfSound.^3) ...
                    .*Prop.Viscosity .*(4/3+(Prop.Gamma-1)./Prop.Prandtl);
        alpha_walls =  2.^(-3/2).*sqrt(omega.*Prop.Viscosity./(Prop.Density.*Prop.SpeedOfSound.^2))...
                      .*(1 + (Prop.Gamma - 1)./sqrt(Prop.Prandtl)).*...
                      Model.Perim./Model.A;

        k =  omega./Prop.SpeedOfSound + (1 - 1i).*alpha_walls;
               
        assignin('base','alpha_wall_WideDuct', alpha_walls)
        alpha_fluid =0;% FluidLosses(omega,Prop);

        WaveNumber.Upstream = k - 1i*alpha_fluid; %The conjugate is taken to be consistent with the definition used throughout the program.
        WaveNumber.Downstream = k - 1i*alpha_fluid; 
        
    case 'RectangularDuct'
        a = Model.a;
        b = Model.b;
        
        %The wave numbers for rectangular ducts are obtained from 
    case 'FluidLosses'
        %Calculate the fluid losses
        alpha_fluid = FluidLosses(omega,Prop);
        
        %Wave number is composed of the 
        %1 Wall absorption 1st order
        %2 Wall phase change 1st order
        %3 Wall absorption 2nd order
        %4 Fluid absorption

        %The wall losses are given by
        r = Model.r;
        Sh = r.*(omega .* Prop.Density./Prop.Viscosity).^(0.5);
        alpha_wall_2ndOrder = omega./Prop.SpeedOfSound.* (1./Sh.^2 .* (1+(Prop.Gamma-1)./sqrt(Prop.Prandtl)-0.5.*Prop.Gamma.*(Prop.Gamma-1)./sqrt(Prop.Prandtl)));

        k_wall_1stOrder = omega./Prop.SpeedOfSound.* ( (1-1i)./(sqrt(2).*Sh).* (1 + (Prop.Gamma-1)./sqrt(Prop.Prandtl)));
     
        k = omega./Prop.SpeedOfSound + 1*k_wall_1stOrder - 1i*alpha_fluid - 1i*alpha_wall_2ndOrder;
        
        WaveNumber.Upstream = real(k) + 1.0*1i*imag(k);
        WaveNumber.Downstream = real(k) + 1.0*1i*imag(k);
end

%%% If the output is called for of this function it is packaged into the output.
if pars.Results.GetOutput
    assignin('base','WaveNumber',WaveNumber);
end

end

function [alpha_fluid] = FluidLosses(omega,Prop)
%Function to calculate the fluid losses
        %This model is taken from "Attenuation of sound in wide ducts with flow
        %at elevated pressure and temperature", C.Lahiri, K.Knobloch, F.Bake,
        %L.Enghardt. Journal of Sound and Vibration 2014 vol 333 page 3440-3458
        f = omega./(2*pi);
        %Equation 15
        alpha_cl = (omega).^2./ (2.*Prop.Density.*Prop.SpeedOfSound.^3).*Prop.Viscosity .*(4/3+(Prop.Gamma-1)./Prop.Prandtl);

        %Equation 16
        R = 287.05;
        Z_rot =  61.6.*exp(-16.8.*(Prop.AmbTemperature+273.15).^(-1/3));
        alpha_rot = (omega).^2 ./ (2.*Prop.Density.*Prop.SpeedOfSound.^3) .* ...
                    Prop.Gamma.*(Prop.Gamma-1).*R./(1.25.*Prop.SpcfHeatCap).*Prop.Viscosity.*Z_rot;

        %Equation 17
        p_ref = 101.325e3;
        T_ref = 293.15;

        T_c = 647.096;
        T_rel = 1-(Prop.AmbTemperature+273.15)./T_c;

        p_ws = 22.064e6.*exp( T_c./(Prop.AmbTemperature+273.15).*(...
            -7.85951783.*T_rel     + ...
             1.84408529.*T_rel.^1.5 + ...
            -11.7866497.*T_rel.^3   + ...
             22.6807411.*T_rel.^3.5 + ...
            -15.9618719.*T_rel.^4   + ... 
             1.80122501.*T_rel.^7.5)    );

        h = Prop.AmbHumidity.*100.*p_ws./Prop.AmbPressure;
        %Equation 19,20 for O2 
        Theta_O = 2239.1;
        X_O = 0.209;
        alpha_lambda_max_O = 2.*pi./35.*X_O.*(Theta_O./(Prop.AmbTemperature+273.15)).^2.*exp(-Theta_O./(Prop.AmbTemperature+273.15));
        f_O = Prop.AmbPressure./p_ref .* (24+4.04.*10.^4.*h.*(0.02+h)./(0.391+h));
        alpha_vib_O = 2.*f./Prop.SpeedOfSound .* alpha_lambda_max_O .* f./f_O ./(1 + (f./f_O).^2);

        %Equation 19,21 for N2
        Theta_N = 3352.9;
        X_N = 0.781;
        alpha_lambda_max_N = 2.*pi./35.*X_N.*(Theta_N./(Prop.AmbTemperature+273.15)).^2.*exp(-Theta_N./(Prop.AmbTemperature+273.15));
        f_N = Prop.AmbPressure./p_ref .* ((Prop.AmbTemperature+273.15)./T_ref).^(-0.5).*(9+280.*h.*exp(-4.170.*( ((Prop.AmbTemperature+273.15)./T_ref).^(-1/3) - 1 )));
        alpha_vib_N = 2.*f./Prop.SpeedOfSound .* alpha_lambda_max_N .* f./f_N ./(1 + (f./f_N).^2);

        alpha_vib = alpha_vib_O + alpha_vib_N;
        alpha_fluid = alpha_cl + alpha_rot+alpha_vib;
end

function [U] = PitotPressureToVelocity(Pressure,Correction,Prop)
%Function to calculate the flow speed based on the measuremt pitot pressure
%and a correction factor based on the flow profile and experimental data
    U = sqrt(2*Pressure./Prop.Density).*Correction;
end
    
