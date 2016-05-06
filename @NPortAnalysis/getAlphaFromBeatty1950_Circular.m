function [Gamma_beatty] = getAlphaFromBeatty1950_Circular(m,n,sh,He,Pr,gamma,RootsBesselFunction)
% beatty 1950 and the erratum
%m is the circumferential mode number
%n is the radial mode number



b = (gamma-1)/sqrt(Pr);

alpha_mn = RootsBesselFunction(m+1,n+1);
He_cuton =  RootsBesselFunction(1,n+1)*pi;

if m == 0
    Gamma_beatty = (1-(He_cuton./He).^2).^(-.5)*...
              1/sqrt(2)*1/sh*(1 + b);
else
Gamma_beatty = (1-(He_cuton./He).^2).^(-.5)*...
              1/sqrt(2)*1/sh*(1-(He_cuton/He)^2*(m/(pi*alpha_mn)) + b);
end

end

