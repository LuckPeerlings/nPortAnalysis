function [Gamma_beatty] = getAlphaFromBeatty1950_Rectangular(ny,nz,sh,He,aspectRatio,Pr,gamma)
% beatty 1950 and the erratum

AR = aspectRatio;
b = (gamma-1)/sqrt(Pr);
He_cuton_y = ny*pi/AR;
He_cuton_z = nz*pi;
He_cuton = sqrt(He_cuton_y^2+He_cuton_z^2);
sin2theta_y = 1-(He_cuton_y./He).^2;
sin2theta_z = 1-(He_cuton_z./He).^2;
coef_y = 2;
coef_z = 2;
if ny==0
    coef_y=1;
end
if nz==0
    coef_z=1;
end
Gamma_beatty = (1-(He_cuton./He).^2).^(-.5)/sqrt(2)./sh.*...
    ((sin2theta_y+b)/AR.*coef_y+(sin2theta_z+b).*coef_z);
end

