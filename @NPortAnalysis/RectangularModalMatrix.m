function ModalMatrix = RectangularModalMatrix(GasProp,WaveNumberProp,Index)

%First calculate which wave numbers will be cut on for the frequency

%The modal matrix is a combination of the ModeShape and the axial
%wavenumber. The Wavenumber in upstream and downstream direction is given
%by k_plus(M) = k(M), k_min(M) = k(-M);

%Model.x is related to the length of the waveguide
%Model.y is related to W the width of the duct and mode m
%Model.z is related to H the height of the duct and mode n
if ~iscolumn(WaveNumberProp.Model.x)
    error('x is not a column vector')
end
if ~iscolumn(WaveNumberProp.Model.y)
    error('y is not a column vector');
end
if ~iscolumn(WaveNumberProp.Model.z)
    error('z is not a column vector');
end

Prop = AirProperties(GasProp);
M = WaveNumberProp.U(Index)/Prop.SpeedOfSound(Index);
omega = 2*pi*WaveNumberProp.f(Index);
%Determine which modes are cut on
ModeVector = MaxMode(   WaveNumberProp.Model.W,...
                        WaveNumberProp.Model.H,...
                        M,...
                        omega,...
                        Prop,...
                        Index);

Phi = ModeShape(    WaveNumberProp.Model.y + WaveNumberProp.Model.W/2,...
                    WaveNumberProp.Model.z + WaveNumberProp.Model.H/2,...
                    WaveNumberProp.Model.W,...
                    WaveNumberProp.Model.H,...
                    M,...
                    ModeVector);
        
[ExpM_plus, ExpM_min] = ExpMatrix(   WaveNumberProp.Model.x,...
                    WaveNumberProp.Model.W,...
                    WaveNumberProp.Model.H,...
                    M,...
                    omega,...
                    ModeVector,...
                    Prop,...
                    Index);

T_plus = Phi.*ExpM_plus;
T_min = Phi.*ExpM_min;

ModalMatrix = [T_plus T_min];
end

function [ExpMatrix_plus, ExpMatrix_min] = ExpMatrix(x,W,H,M,omega,ModeVector,GasProp,Index)
    k0 = omega/GasProp.SpeedOfSound(Index);
    if M ~= 0
        error('The functions is not yet defined for flows unequal to zero')
    end
    for ii = 1:size(ModeVector,1)
        m = ModeVector(ii,1);
        n = ModeVector(ii,2);
        shV = H*sqrt(GasProp.Density(Index)*omega/GasProp.Viscosity(Index));
        [Gamma_beatty] = NPortAnalysis.getAlphaFromBeatty1950_Rectangular(m,n,shV,...
            k0*H,W/H,GasProp.Prandtl(Index),GasProp.Gamma(Index));
        kx_plus = sqrt_complex(k0.^2 - (m*pi/W)^2 - (n*pi/H)^2 ) + (1-1i)*Gamma_beatty*k0;
        kx_min = -kx_plus;
        ExpMatrix_plus(:,ii) = exp(-1i.*kx_plus.*x);
        ExpMatrix_min(:,ii) = exp(-1i.*kx_min.*x);
    end    
end


function ModeVector = MaxMode(W,H,M,omega,GasProp,Index)
    MaxM = 0;
    MaxN = 0;
    k0 = omega/GasProp.SpeedOfSound(Index);    
    ii = 1;
    for mm = 0:MaxM;
        for nn = 0:MaxN;           
            k_2 = k0^2 - (mm*pi/W)^2 - (nn*pi/H)^2;
            if k_2 >= -0.05*k0^2
                ModeVector(ii,1) = mm;
                ModeVector(ii,2) = nn;
                ii = ii+1;
            end
        end
    end
    assignin('base','ModeVector',ModeVector)
end  
    
function Phi = ModeShape(y,z,W,H,M,ModeVector)
    %The positions should be row vectors,
    assignin('base','y',y)
    assignin('base','z',z)
    for ii = 1:size(ModeVector,1)
        m = ModeVector(ii,1);
        n = ModeVector(ii,2);
        if M > 0
            Phi(:,ii) = cos(m.*pi./W.*y).*cos(n.*pi./H.*z);
        elseif M < 0
            Phi(:,ii) = cos(m.*pi./W.*y).*cos(n.*pi./H.*z);
        else
            Phi(:,ii) = cos(m.*pi./W.*y).*cos(n.*pi./H.*z);
        end        
    end
    assignin('base','Phi',Phi)
end

function Out = sqrt_complex(In)
%This function calculates the square root of a complex number such that the
%imaginary part is always negative.
Out = -sign(imag(sqrt(In)))*sqrt(In);
if imag(sqrt(In)) == 0
    Out = sqrt(In);
end
end