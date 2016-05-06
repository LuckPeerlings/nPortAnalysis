function ModalMatrix = CircularModalMatrix(GasProp,WaveNumberProp,Index)

RootsBesselFunction =   [   0                   3.83170597020751 	7.01558666981561 	10.1734681350627 	13.3236919363142;
                            1.84118378134065 	5.33144277352503 	8.53631636634628 	11.7060049025920 	14.8635886339090;
                            3.05423692822714 	6.70613319415845 	9.96946782308759 	13.1703708560161 	16.3475223183217;
                            4.20118894121052 	8.01523659837595 	11.3459243107430 	14.5858482861670 	17.7887478660664;
                            5.31755312608399 	9.28239628524161 	12.6819084426388 	15.9641070377315 	19.1960288000489;
                            6.41561637570024 	10.5198608737723 	13.9871886301403 	17.3128424878846 	20.5755145213868;
                            7.50126614468414 	11.7349359530427 	15.2681814610978 	18.6374430096662 	21.9317150178022;
                            8.57783648971407 	12.9323862370895 	16.5293658843669 	19.9418533665273 	23.2680529264575;
                            9.64742165199721 	14.1155189078946 	17.7740123669152 	21.2290626228531 	24.5871974863176;
                            10.7114339706999 	15.2867376673329 	19.0045935379460 	22.5013987267772 	25.8912772768391;
                            11.7708766749555 	16.4478527484865 	20.2230314126817 	23.7607158603274 	27.1820215271905;
                          ];
%First calculate which wave numbers will be cut on for the frequency

%The modal matrix is a combination of the ModeShape and the axial
%wavenumber. The Wavenumber in upstream and downstream direction is given
%by k_plus(M) = k(M), k_min(M) = k(-M);

%Model.x is related to the length of the waveguide
%Model.y is related to W the width of the duct and mode m
%Model.z is related to H the width of the duct and mode n
if ~iscolumn(WaveNumberProp.Model.x)
    error('x is not a column vector')
end
if ~iscolumn(WaveNumberProp.Model.r)
    error('y is not a column vector');
end
if ~iscolumn(WaveNumberProp.Model.theta)
    error('theta is not a column vector');
end

Prop = AirProperties(GasProp);
M = WaveNumberProp.U(Index)/Prop.SpeedOfSound(Index);
omega = 2*pi*WaveNumberProp.f(Index);
%Determine which modes are cut on
ModeVector = MaxMode(   WaveNumberProp.Model.R,...
                        M,...
                        omega,...
                        Prop,...
                        RootsBesselFunction,...
                        Index);
                    
Phi = ModeShape(    WaveNumberProp.Model.r,...
                    WaveNumberProp.Model.theta,...
                    WaveNumberProp.Model.R,...
                    M,...
                    ModeVector,...
                    RootsBesselFunction);
        
ExpM = ExpMatrix(   WaveNumberProp.Model.x,...
                    WaveNumberProp.Model.R,...
                    M,...
                    omega,...
                    ModeVector,...
                    Prop,...
                    RootsBesselFunction,...
                    Index);
T_plus = Phi.*ExpM;

Phi = ModeShape(    WaveNumberProp.Model.r,...
                    WaveNumberProp.Model.theta,...
                    WaveNumberProp.Model.R,...
                    -M,...
                    ModeVector,...
                    RootsBesselFunction);
        
ExpM = ExpMatrix(   WaveNumberProp.Model.x,...
                    WaveNumberProp.Model.R,...
                    -M,...
                    omega,...
                    ModeVector,...
                    Prop,...
                    RootsBesselFunction,...
                    Index);
                
T_min = Phi./ExpM;

ModalMatrix = [T_plus T_min];
end

function ExpMatrix = ExpMatrix(x,R,M,omega,ModeVector,GasProp,RootsBesselFunction,Index)
    k0 = omega/GasProp.SpeedOfSound(Index);
    for ii = 1:1;%size(ModeVector,1)
        m = ModeVector(ii,1);
        n = ModeVector(ii,2);
        k_mn = RootsBesselFunction(m+1,n+1)/R;
        
        shV = R*sqrt(GasProp.Density(Index)*omega/GasProp.Viscosity(Index));
        [Gamma_beatty] = NPortAnalysis.getAlphaFromBeatty1950_Circular(m,n,shV,...
            k0*R,GasProp.Prandtl(Index),GasProp.Gamma(Index),RootsBesselFunction);        
       
        kx = (k0*M + sqrt(k0^2 - (1-M^2)*k_mn^2))/(1-M^2);             
        k_tot = kx + (1-1i)*Gamma_beatty*k0/(1+M);           
        if M == 0
        	ExpMatrix(:,ii) = exp(-1i.*k_tot.*x);
            
        else
            ExpMatrix(:,ii) = exp(sign(M)*1i.*k_tot.*x);
        end
    end    
end


function ModeVector = MaxMode(R,M,omega,GasProp,RootsBesselFunction,Index)
    MaxM = 3;
    MaxN = 3;
    ii = 1;
    k0 = omega/GasProp.SpeedOfSound(Index);
    for mm = 0:MaxM;
        for nn = 0:MaxN;
            k_mn = RootsBesselFunction(mm+1,nn+1)/R;
            k_2 = k0^2 - (1-M^2)*(k_mn)^2;          
            if k_2 > 0
                ModeVector(ii,1) = mm;
                ModeVector(ii,2) = nn;
                ii = ii+1;
            end
        end
    end
end  
    
function Phi = ModeShape(r,theta,R,M,ModeVector,RootsBesselFunction)
    %The positions should be row vectors,
    for ii = 1:size(ModeVector,1)
        m = ModeVector(ii,1);
        n = ModeVector(ii,2);
        if M > 0
            Phi(:,ii) = besselj(m,r/R.*RootsBesselFunction(m+1,n+1)).*exp(1i.*theta.*m);
        elseif M < 0
            Phi(:,ii) = besselj(m,r/R.*RootsBesselFunction(m+1,n+1)).*exp(1i.*theta.*m);
        else
            Phi(:,ii) = besselj(m,r/R.*RootsBesselFunction(m+1,n+1)).*exp(1i.*theta.*m);
        end        
    end
end
