function [x_avg] = MicPosCalibration(f,Port1,fig)

if isfield(Port1,'Constant')
   Port1.Meas1 = MergeStructs(Port1.Meas1,Port1.Constant);
end
                    
Z = Port1.Meas1.P;
k = WaveNumber(f,Port1);
x_mic = Port1.Meas1.x;
%Calculating the transfer functions between the microphones for each two
%pair of microphones
COUNTER = zeros(1,length(x_mic));
XPos = zeros(1,length(x_mic));


ss = 1;
%Determine the best fit for each microphone combination
for ii = 1:length(x_mic)
    for jj =ii+1:length(x_mic)
    x0 = [x_mic(ii),x_mic(jj)];
    x0s(:,ss) = x0;
    DATA = Z(jj,:)./Z(ii,:);
    options = optimset('Display', 'off','LargeScale','off') ;
    x_calc = fminunc(@(x) MicTF(x,DATA,k,f,Port1),x0,options);
    x_calc_s(:,ss) = x_calc;
    ss = ss+1;
    XPos(ii) = XPos(ii) + x_calc(1);
    COUNTER(ii) = COUNTER(ii)+1;
    XPos(jj) = XPos(jj) + x_calc(2);
    COUNTER(jj) = COUNTER(jj)+1;
    end
end
%Determine the "average" best fit
x_avg = XPos./COUNTER;
if fig == 1
figure
hold all
ylabel(gca, 'Re ( P(2)/P(1))',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
xlabel(gca, 'Frequency',...
       'Fontsize',20, ...
       'Fontweight', 'bold'...
       );

set(gca, ...
        'Box'           , 'off'     , ...
        'TickDir'       , 'out'     , ...                
        'XScale'        ,'linear'      ,...      
        'Fontsize'      ,16         ,...
        'FontWeight'    ,'Bold',...
        'FontName','Times',...
        'XMinorTick'    ,'on'       ,...
        'YMinorTick'    ,'on'       ,...
        'XColor'      , [0 0 0], ...
        'YColor'      , [0 0 0], ...
        'LineWidth'     ,0.5,...
        'Xgrid'         ,'off',...
        'Ygrid'         ,'off',... 
        'Box'           ,'on'...
    );
set(gcf, 'color', 'white');

R = ReflCoeff(f,Port1);
TF_0   = ( exp(-1i.*k.*x_mic(2)) + R.*exp(1i.*k.*x_mic(2)))./( exp(-1i.*k.*x_mic(1)) + R.*exp(1i.*k.*x_mic(1)));
TF_opt = ( exp(-1i.*k.*x_avg(2)) + R.*exp(1i.*k.*x_avg(2)))./( exp(-1i.*k.*x_avg(1)) + R.*exp(1i.*k.*x_avg(1)));
TF_Meas = Z(2,:)./Z(1,:);
plot(f,real(TF_Meas),'k')
plot(f,real(TF_opt))
plot(f,real(TF_0))
legend('Measured','Optimized Solution','Starting Value')
end
end

function err = MicTF(x,DATA,k,f,Port1)
R = ReflCoeff(f,Port1);
model = ( exp(-1i.*k.*x(2)) + R.*exp(1i.*k.*x(2)))./( exp(-1i.*k.*x(1)) + R.*exp(1i.*k.*x(1)));
err = sum(abs(DATA./model - 1));
end

function R = ReflCoeff(f,Port1)
    Prop = GasProperties(Port1.Meas1.GasProp);
    omega = 2.*pi.*f;
    ZT = 1./(exp(1i*pi/4)./(Prop.Density.*Prop.SpeedOfSound).*sqrt(omega.*Prop.Viscosity./(Prop.Density.*Prop.SpeedOfSound.^2)).*(Prop.Gamma-1)./sqrt(Prop.Prandtl));

    Zs = Prop.Density .* Prop.SpeedOfSound;
    R = (ZT - Zs) ./ (Zs+ZT);
end

function k = WaveNumber(f,Port1)
    %Calculate the fluid losses
    Prop = GasProperties(Port1.Meas1.GasProp);
    r = Port1.Meas1.WaveNumberProp.Model.r;
    alpha_fluid = FluidLosses(f,Prop);
    omega = 2.*pi.*f;
    
    %Wave number is composed of the 
    %1 Wall absorption 1st order
    %2 Wall phase change 1st order
    %3 Wall absorption 2nd order
    %4 Fluid absorption

    %The wall losses are given by
    Sh = r.*(omega .* Prop.Density./Prop.Viscosity).^(0.5);
    alpha_wall_2ndOrder = omega./Prop.SpeedOfSound.* (1./Sh.^2 .* (1+(Prop.Gamma-1)./sqrt(Prop.Prandtl)-0.5.*Prop.Gamma.*(Prop.Gamma-1)./sqrt(Prop.Prandtl)));

    k_wall_1stOrder = omega./Prop.SpeedOfSound.* ( (1-1i)./(sqrt(2).*Sh).* (1 + (Prop.Gamma-1)./sqrt(Prop.Prandtl)));

    k = omega./Prop.SpeedOfSound + 1*k_wall_1stOrder - 1i*alpha_fluid - 1i*alpha_wall_2ndOrder;
end  

function [alpha_fluid] = FluidLosses(f,Prop)
%Function to calculate the fluid losses
        %This model is taken from "Attenuation of sound in wide ducts with flow
        %at elevated pressure and temperature", C.Lahiri, K.Knobloch, F.Bake,
        %L.Enghardt. Journal of Sound and Vibration 2014 vol 333 page 3440-3458
        omega = 2*pi*f;
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
    

