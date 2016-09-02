classdef AreaExpansionModel  < matlab.mixin.SetGet

    properties 
        FreqVec
        RadiusUpstream
        RadiusDownstream
        Temperature
        Humidity
        AmbientPressure
        ModelName
        ScatMatrix
    end
    methods
        function obj = CalculateScatteringMatrix(obj)
        
        AirProp.t = obj.Temperature;
        AirProp.RH = obj.Humidity;
        AirProp.p = obj.AmbientPressure;
        
        Prop = AirProperties(AirProp);
        switch obj.ModelName
            case 'Kergomard'
            Zj = GetImpedanceKergoMard(obj,2*pi*obj.FreqVec,Prop);
            case 'Auregan'
            Zj = GetImpedanceAuregain(obj,2*pi*obj.FreqVec,Prop);                   
        end
        
        mu = (obj.RadiusUpstream/obj.RadiusDownstream).^2;
        Sa = pi*(obj.RadiusUpstream)^2;
        
        % Using the equation 7-3.5 pag 330 of "Acoustics: An Introduction to Its
        % Physical Principels and Applications, Allan D. Pierce
        obj.ScatMatrix.S11 = 1./(mu + 1 + Sa/(Prop.Density*Prop.SpeedOfSound)*Zj) .* (mu - 1 + Sa/(Prop.Density*Prop.SpeedOfSound)*Zj);
        obj.ScatMatrix.S12 = 1./(mu + 1 + Sa/(Prop.Density*Prop.SpeedOfSound)*Zj ) .* 2;
        obj.ScatMatrix.S21 = 1./(mu + 1 + Sa/(Prop.Density*Prop.SpeedOfSound)*Zj ) .* (2.*mu );
        obj.ScatMatrix.S22 = 1./(mu + 1 + Sa/(Prop.Density*Prop.SpeedOfSound)*Zj ) .* (1 - mu + Sa/(Prop.Density*Prop.SpeedOfSound)*Zj);
        end
        
        function Zj = GetImpedanceKergoMard(obj,omega,Prop)
            alpha = obj.RadiusUpstream/obj.RadiusDownstream;
            eta = 1-alpha;
            if alpha < 0.55
                L0 =  0.26164 - 0.353*alpha + 0.0809*alpha.^3 + 0.00119*alpha.^5 + 0.0175*alpha.^6;
                A = 1.8598.*eta.^3 - 1.2029.*eta.^4 + 0.162.*eta.^5;
                B = -0.1197 + 0.6488.*alpha - 1.074.*alpha.^2 + 0.85.*alpha.^4;
            else
                L0 = Prop.Density./obj.RadiusUpstream.*(4/pi^2).*eta.^2.*(-0.49198.*log(eta) + 0.50349 - 0.376246*eta.^2 - 0.852222.*eta.^2.*log(eta));
                A = 0.0000137 + 0.1889.*eta.^2 + 0.8202.*eta.^3 + 0.9301.*eta.^4 - 1.47586.*eta.^5;
                B = -0.01554.*eta.^2 - 0.4249.*eta.^2 + 0.776.*eta.^4;
            end
            k_c = 3.831706./obj.RadiusDownstream;
            Q = 1./(sqrt(1-((omega./Prop.SpeedOfSound)./k_c).^2)) - 1;
            
            Zj = (1i.*omega).*(L0 + alpha.*Prop.Density./obj.RadiusUpstream.*A.*Q + alpha.*Prop.Density./obj.RadiusUpstream.*B.*Q.^2) ;                     
        end
        
        function Zj = GetImpedanceAuregain(obj,omega,Prop)
            Sa = pi*obj.RadiusUpstream^2 ;
            alpha = obj.RadiusUpstream/obj.RadiusDownstream;
          
            f_alpha = ( (1-alpha)*(3+alpha) )/( 3*(1+alpha)^2);
            gamma = 2 / (alpha^2*(1-alpha^2)*(0.25+f_alpha));
            k2 = ((omega/Prop.SpeedOfSound).^2 - (gamma/obj.RadiusDownstream.^2)).^(0.5);
            dL = (1-alpha.^2)./(1i.*k2);
            Zj = (-1i.*omega .* Prop.Density/Sa).*dL;
        end
        
        function displayScatMatrix(obj)
            figure;            
            for ii = 1:2
                for jj = 1:2
                    subplot(2,2, (ii - 1)*2 + jj )
                    plot(obj.FreqVec, abs( obj.ScatMatrix.(['S',num2str(jj),num2str(ii)]) ))
                end
            end
            figure;            
            for ii = 1:2
                for jj = 1:2
                    subplot(2,2, (ii - 1)*2 + jj )
                    plot(obj.FreqVec, unwrap(angle( obj.ScatMatrix.(['S',num2str(jj),num2str(ii)]) ))*180/pi)
                end
            end
        end
    end
end