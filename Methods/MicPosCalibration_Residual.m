function [x_avg, T_corr] = MicPosCalibration_Residual(f,Port1,MeasSelect,fig)

if isfield(Port1,'Constant')
   Port1.Meas1 = MergeStructs(Port1.Meas1,Port1.Constant);
end
Z = Port1.Meas1.P(:,MeasSelect);
f = f(MeasSelect);
Port1.Meas1.GasProp.t = Port1.Meas1.GasProp.t(MeasSelect);
x_mic = Port1.Meas1.x;


x0 = [x_mic,0];
ResTot = GetResidual(x0,f,Port1,Z,false,MeasSelect);

fun = @(x) GetResidual(x,f,Port1,Z,false,MeasSelect);
lb(1:length(x_mic))= x_mic-10e-3;
ub(1:length(x_mic))= x_mic+10e-3;
lb(end+1) = 0;
ub(end+1) = 2;


options = optimoptions(@fminunc,'Display','iter-detailed','MaxFunEvals',1000,'TolFun',1e-15,'TolX',1e-20,'Maxiter',10000);
x = fmincon(fun,x0,[],[],[],[],lb,ub,[]);
%x1 = fminunc(fun,x0,options)
x0;
x_avg = x(1:end-1);
T_corr = x(end);



fun = @(x) GetResidual(x,f,Port1,Z,true,MeasSelect);
fun(x0)
fun(x)

end

function ResTot = GetResidual(x,f,Port1,Z,ShowSingleRes,MeasSelect)
%k = k0*1/x(end);
Port1.Meas1.GasProp.t = Port1.Meas1.GasProp.t+x(end);

k = WaveNumber(f,Port1);
A = GetModalMatrix(x(1:end-1),k);
%Calculate the combined residual
for ff = 1:length(k)
    A_single = squeeze(A(:,:,ff));
    P_decomp = pinv(A_single)*Z(:,ff);
    P_model = A_single*P_decomp;
    res(ff) = (P_model-Z(:,ff))'*(P_model-Z(:,ff));
end
if ShowSingleRes
    figure; plot(MeasSelect,res);
end

ResTot = sum(res) + sum(abs(imag(P_decomp(1,:)./P_decomp(2,:))));
end

function A = GetModalMatrix(x_mic,k)
for ff = 1:length(k)
    for ll = 1:length(x_mic);
        A(ll,1,ff) = exp(-1i*k(ff)*x_mic(ll));
        A(ll,2,ff) = exp(+1i*k(ff)*x_mic(ll));
    end
end
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
    

