%TODO List

%Create function that calculates the k_y for a given frequency and
%admittance

%Determine the resulting waves for a given impedance, using the
%conservation of mass and momentum and create function to prove it.


classdef ImpedanceEduction 
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ad
        
    end
    properties (SetAccess = private, GetAccess = public)
        %FluidProperties %Variable containing the properties of the fluid in the duct
        %Omega
        %DuctHeigth
    end
    
    methods
%         function ImpedanceEduction(obj)
%             
%         end
        
        function TestTransversalWaveNumber(obj)    
            %Test function to determine the proper working if the
            %determination of the axial wave number. 
            
            %The hard wall case is checked, for which the analytical
            %solutions are given by k_y = n*pi/(2*DuctHeight)
            
            %Depending on the initial guess, different solutions will be
            %found, corresponding to the cases n = 1,2,3,...
            
            %Initial values
            Y_ac_Init = 1+1i;
            M = 1e-4;
            k_x = 0;
            Omega_vec = linspace(1e-7,10,100)*2*pi;
            FluidProperties.SpeedOfSound = 343;
            FluidProperties.Density = 1.7;
            DuctHeight = 15e-3;
            NrOfSteps = 1000;
            
            %Using two initial guesses
            InitialGuess_ky_0 = 0*pi/(2*DuctHeight);
            InitialGuess_ky_1 = 1*pi/(2*DuctHeight);
            InitialGuess_ky_2 = 2*pi/(2*DuctHeight);
            
            %Solving the transversal wave numbers
            for ff = 1:length(Omega_vec)
                %Only
                if ff ~= 1
                    NrOfSteps = 1;
                end
                k_y_0(ff) = ImpedanceEduction.SolveTranversalWaveNumber(Y_ac_Init, InitialGuess_ky_0, M, k_x, NrOfSteps, DuctHeight,Omega_vec(ff),FluidProperties);
                k_y_1(ff) = ImpedanceEduction.SolveTranversalWaveNumber(Y_ac_Init, InitialGuess_ky_1, M, k_x, NrOfSteps, DuctHeight,Omega_vec(ff),FluidProperties);
                k_y_2(ff) = ImpedanceEduction.SolveTranversalWaveNumber(Y_ac_Init, InitialGuess_ky_2, M, k_x, NrOfSteps, DuctHeight,Omega_vec(ff),FluidProperties);
            end            
            %Plotting the results
            figure;
            subplot(2,1,1); 
                plot(Omega_vec/(2*pi),real(k_y_0),'x');
                hold all
                plot(Omega_vec/(2*pi),real(k_y_1),'x');
                plot(Omega_vec/(2*pi),real(k_y_2),'x');
                n = 0; plot(Omega_vec/(2*pi),n*pi/(2*DuctHeight)*ones(size(Omega_vec)));
                n = 1; plot(Omega_vec/(2*pi),n*pi/(2*DuctHeight)*ones(size(Omega_vec)));
                n = 2; plot(Omega_vec/(2*pi),n*pi/(2*DuctHeight)*ones(size(Omega_vec)));
                legend('Optimized value: Starting Value 0','Optimized value: Starting Value 1','Optimized value: Starting Value 2','Mode 0','Mode 1','Mode 2')
                xlabel('Frequency [Hz]')
                ylabel('Real part of the wave number')
                title('Transversal wavenumbers')
            subplot(2,1,2);
                plot(Omega_vec/(2*pi),imag(k_y_0),'x');
                hold all
                plot(Omega_vec/(2*pi),imag(k_y_1),'x');
                plot(Omega_vec/(2*pi),imag(k_y_2),'x');
                plot(Omega_vec/(2*pi),0*ones(size(Omega_vec)));
                xlabel('Frequency [Hz]')
                ylabel('Imaginary part of the wave number')
        end
        
        function TestModeMatching(obj)
            %This function is made to check the Mode matching method. The
            %modes in a hard walled duct are calculated for a hard walled
            %section
            P_plus_a = 1;
            P_plus_c = 0.02;
            NrModes = 5;
            NrOfSteps_WaveNumberCalculation = 100;
            Y_ac_Init = 1e2 + 1i*1e2;
            H = 1;
            M = 0.1;
            Omega = 1;
            FluidProperties.Density = 1.7;
            FluidProperties.SpeedOfSound = 343*(1+1e-2*1i);
            LinerLength=0.50;
            Amplitudes = ImpedanceEduction.ModeMatching(P_plus_a,P_plus_c,H,LinerLength,Omega,M,FluidProperties,NrModes,NrOfSteps_WaveNumberCalculation,Y_ac_Init);
            fprintf('----------------------------------------------------\n')
            fprintf('Hard walled test case, mean flow speed of M %1.3f \n',M)
            fprintf('Inputs \t \t \t \t \t \t Outputs \n')
            fprintf('P_Plus_a %6.2f+1i*%6.2f \t P_min_c %6.2f+1i*%6.2f \n', real(P_plus_a),imag(P_plus_a)  , real(Amplitudes(2)), imag(Amplitudes(2)))
            fprintf('P_Plus_c %6.2f+1i*%6.2f \t P_min_a %6.2f+1i*%6.2f \n', real(P_plus_c),imag(P_plus_c) , real(Amplitudes(1)), imag(Amplitudes(1)))                       
            fprintf('\nIntermediate waves \n')
            for nn = 1:NrModes
                fprintf('Mode %i \t P_Plus_b %5.2f+1i*%5.2f \t  P_Min_b %5.2f+1i*%5.2f \n',nn, ...
                            real(Amplitudes((nn-1)*4 + 3 )), imag(Amplitudes((nn-1)*4 + 3)),...
                            real(Amplitudes((nn-1)*4 + 4)), imag(Amplitudes((nn-1)*4 + 4)))
                              
            end
        end
    end
    
    methods(Static)
        function [Amplitudes] = ModeMatching(P_plus_a,P_plus_c,DuctHeight,LinerLength,Omega,MachNumber,FluidProperties,NrModes,NrOfSteps_WaveNumberCalculation,Y_ac_Init)
            if LinerLength <= 0
                warning('Negative or zero liner length, setting liner length to 1e-6m');
                LinerLength = 1e-6;
            end
            %Returns the amplitudes from the ModeMatching for the N modes. 
            %The results are ordered by
            %Amplitudes = [P_min_a, P_min_c, P_plus_b, P_min_b]
            %and each mode is concatonated, starting with mode 0. 
            %Schematic of the wave definitions
            %--------------------------------------------------------------
            %    --> P_plus_a  |       ---> P_plus_b  |   <--- P_plus_c
            %                  |                      |
            %                  |          ------\     |
            %                  | Flow Direction  >    |
            %                  |          ------/     |
            %                  |                      |
            %    <-- P_min_a   |      <--- P_min_b    |   ---> P_min_c
            %--------------------------------------------------------------
            %                x=0                      x=L
            %
            %   /\ k_y
            %    |
            %k_x 0--> k_z
            %
            %This function performs the mode matching to obtain the
            %pressure amplitudes of the waves, based on the given
            %information
            
            %For these calculations, the wavenumber in the third direction
            %is assumed to be only corresponding to the plane wave mode
            K_X_Plane = 0;   
            k0 = Omega/FluidProperties.SpeedOfSound; %Free field wave number
          
            %Calculate the transveral wave numbers of the N-modes in the HardWalled section;
            for nn = 1:NrModes
                k_y_HardWalledSec(nn) = pi*(nn-1)/(DuctHeight);
            end
            %Calculate the transveral wave numbers of the M-modes in the lined section; 
            for mm = 1:NrModes
                k_y_LinedSec(mm) = ImpedanceEduction.SolveTranversalWaveNumber(Y_ac_Init, (mm-1)*pi/(DuctHeight), MachNumber, K_X_Plane, NrOfSteps_WaveNumberCalculation, DuctHeight,Omega,FluidProperties);
            end
            
            %Calculate the axial wave-numbers in the hardwalled section and
            %the lined section
            for nn = 1:NrModes
                [k_axial_HardWalledSec.WithFlow(nn), k_axial_HardWalledSec.AgainstFlow(nn)] = ...
                    ImpedanceEduction.AxialWaveNumber( MachNumber,k0,K_X_Plane,k_y_HardWalledSec(nn));
                
                [k_axial_LinedSec.WithFlow(nn), k_axial_LinedSec.AgainstFlow(nn)] = ...
                    ImpedanceEduction.AxialWaveNumber( MachNumber,k0,K_X_Plane,k_y_LinedSec(nn));
            end
            
            %Calculate the factors that arise when the complete solution is
            %multiplied by the hard walled mode shapes and integrated over
            %the duct height.
            for nn = 1:NrModes
                for mm = 1:NrModes
                    %Calculation of delta for the first duct (hard-walled)
                    if nn == mm
                        Delta_11(nn,mm) = DuctHeight/2;
                        if nn == 1 && mm == 1
                            Delta_11(nn,mm) = DuctHeight;
                        end
                    end
                    %Calculation of delta for the third duct (hard-walled)
                    if nn == mm
                        Delta_13(nn,mm) = DuctHeight/2;
                        if nn == 1 && mm == 1
                            Delta_13(nn,mm) = DuctHeight;
                        end
                    end
                    %Calculation of delta for the second duct (lined)
                    Delta_12(nn,mm) = 1/2*( sin( (k_y_HardWalledSec(nn) + k_y_LinedSec(mm) )*DuctHeight)/(k_y_HardWalledSec(nn) + k_y_LinedSec(mm)) ...
                        + sin((k_y_HardWalledSec(nn) - k_y_LinedSec(mm))*DuctHeight)/(k_y_HardWalledSec(nn) - k_y_LinedSec(mm)));
                    %When the transversal wave number of the hard
                    %walled and lined section are equal,  the
                    %denominator of the second term becomes zero. Next
                    %part is to handle that exception.
                    if k_y_HardWalledSec(nn) == k_y_LinedSec(mm)
                        Delta_12(nn,mm) = 1/2*( sin((k_y_HardWalledSec(nn) + k_y_LinedSec(mm))*DuctHeight)/(k_y_HardWalledSec(nn) + k_y_LinedSec(mm)) ...
                            + DuctHeight);
                        if k_y_HardWalledSec(nn) == 0
                            Delta_12(nn,mm) = DuctHeight;
                        end
                    end
                 end
            end
            %The solution vector contains is expressed as
            %[P_min_a, P_min_c, P_plus_b, P_min_b]
            %Corresponding to the Theory as
            %[B^1 A^3 A^2 B^2]
            %For each mode
            %Create the matrices for the mode matching
            
            %For each mode, there are four equations
            %Only for the first mode, these equations are inhomogeneous
            A = zeros(4*NrModes);
            for nn = 1:NrModes
                for mm = 1:NrModes
                    %Boundary Condition (1) Continuity of pressure at x=0
                    A((nn-1)*4 + 1,(mm-1)*4 + 1) = -Delta_11(nn,mm);
                    A((nn-1)*4 + 1,(mm-1)*4 + 3) =  Delta_12(nn,mm);
                    A((nn-1)*4 + 1,(mm-1)*4 + 4) =  Delta_12(nn,mm);
                    
                    %Boundary Condition (2) Continuity of velocity at x=0
                    A((nn-1)*4 + 2,(mm-1)*4 + 1) = -Delta_11(nn,mm)*k_axial_HardWalledSec.AgainstFlow(mm);
                    A((nn-1)*4 + 2,(mm-1)*4 + 3) =  Delta_12(nn,mm)*k_axial_LinedSec.WithFlow(mm);
                    A((nn-1)*4 + 2,(mm-1)*4 + 4) =  Delta_12(nn,mm)*k_axial_LinedSec.AgainstFlow(mm);
                    
                    %Boundary Condition (3) Continuity of pressure at x=L
                    A((nn-1)*4 + 3,(mm-1)*4 + 2) = ...
                            -Delta_13(nn,mm)*exp(-1i*k_axial_HardWalledSec.WithFlow(mm)*LinerLength) ;
                    A((nn-1)*4 + 3,(mm-1)*4 + 3) = ...
                             Delta_12(nn,mm)*exp(-1i*k_axial_LinedSec.WithFlow(mm)*LinerLength);
                    A((nn-1)*4 + 3,(mm-1)*4 + 4) = ...
                             Delta_12(nn,mm)*exp(-1i*k_axial_LinedSec.AgainstFlow(mm)*LinerLength);
                    
                    %Boundary Condition (4) Continuity of velocity at x=L
                    A((nn-1)*4 + 4,(mm-1)*4 + 2) = ...
                            -Delta_13(nn,mm)*k_axial_HardWalledSec.WithFlow(mm)*exp(-1i*k_axial_HardWalledSec.WithFlow(mm)*LinerLength) ;
                    A((nn-1)*4 + 4,(mm-1)*4 + 3) = ...
                             Delta_12(nn,mm)*k_axial_LinedSec.WithFlow(mm)*exp(-1i*k_axial_LinedSec.WithFlow(mm)*LinerLength);
                    A((nn-1)*4 + 4,(mm-1)*4 + 4) = ...
                             Delta_12(nn,mm)*k_axial_LinedSec.AgainstFlow(mm)*exp(-1i*k_axial_LinedSec.AgainstFlow(mm)*LinerLength);
                end
            end
            %Create the vector of boundary conditions
            B = zeros(4*NrModes,1);
            %This vector only contains the plane wave mode, mode 1,
            nn = 1;
            for mm = 1:NrModes
                B((mm-1)*4 + 1) = Delta_11(nn,mm)*P_plus_a ;
                B((mm-1)*4 + 2) = Delta_11(nn,mm)*P_plus_a*k_axial_HardWalledSec.WithFlow(mm);
                B((mm-1)*4 + 3) = Delta_13(nn,mm)*P_plus_c*                                        exp(-1i*k_axial_HardWalledSec.AgainstFlow(mm)*LinerLength);
                B((mm-1)*4 + 4) = Delta_13(nn,mm)*P_plus_c*k_axial_HardWalledSec.AgainstFlow(mm)*  exp(-1i*k_axial_HardWalledSec.AgainstFlow(mm)*LinerLength);
            end
            %Solve the system of equations
            Amplitudes = A\B;
            
        end
       
        function Y_ac = SolveImpedance(InitialGuess_Y_ac, ky_Init, M, k_x, NrOfSteps, DuctHeigth,Omega,FluidProperties)
            %This function solves the EigenValueEquation by taking a suitable
            %initial guess of the acoustic impedance and changing the
            %transversal wavenumber and Mach number to determine the acoustic impedance

            %The stepping can have a logarithmic increase
            M_Vec = linspace(0,1,NrOfSteps).^(1/1)*M;
            ky_Vec = linspace(0,1,NrOfSteps).^(1/1)*ky_Init;
            
            %Initialize the wavenumber vector in the y-direction and step
            %towards the final acoustic impedance and Mach Number
            Y_ac_vec(1) = InitialGuess_Y_ac;
            for ii  = 2:NrOfSteps
                func = @(Y_ac) abs(ImpedanceEduction.EigenValueEquation(Y_ac,M_Vec(ii),k_x,ky_Vec(ii),DuctHeigth,Omega,FluidProperties));
                Y_ac_vec(ii) = fsolve(func,Y_ac_vec(ii-1));
            end
            
            Y_ac =  Y_ac_vec(end);                
        end
        
        function k_y = SolveTranversalWaveNumber(Y_ac_Init, InitialGuess_ky, M, k_x, NrOfSteps, DuctHeigth,Omega,FluidProperties)
            %This function solves the EigenValueEquation by taking a suitable
            %initial guess of the transversal wave number and changing the
            %impedance and Mach number to determine the transveral
            %wavenumber
            
            %If the NrOfSteps is equal to 1, there will be no steps and the
            %iteration starts directly from the given flow number and
            %admittance

            %The stepping can have a logarithmic increase
            M_Vec = linspace(0,1,NrOfSteps).^(1/1)*M;
            Y_Vec = linspace(0,1,NrOfSteps).^(1/1)*Y_ac_Init;
            
            %Initialize the wavenumber vector in the y-direction and step
            %towards the final acoustic impedance and Mach Number
            k_y_vec(1) = InitialGuess_ky;
            options = optimoptions('fsolve','FiniteDifferenceType','Central','Display','Iter-detailed');
            for ii  = 1:NrOfSteps
                func = @(k_y) ImpedanceEduction.EigenValueEquation(Y_Vec(ii),M_Vec(ii),k_x,k_y,DuctHeigth,Omega,FluidProperties);
                k_y_vec(ii+1) = fsolve(func,k_y_vec(ii),options);
            end            
            k_y = k_y_vec(end);                
        end
        
        function Res = EigenValueEquation(Y_ac,M,k_x,k_y,DuctHeight,Omega,FluidProperties)
            %The eigenvalue equation that has to be solved to obtain the
            %wavenumbers or acoustic impedance
            k0 = Omega/FluidProperties.SpeedOfSound;
            [AxialWaveNumber_withFlow,~] = ImpedanceEduction.AxialWaveNumber(M,k0,k_x,k_y);
            Res = k_y*tan(2*k_y*DuctHeight) - 1i*FluidProperties.Density*FluidProperties.SpeedOfSound^2*Y_ac/Omega *(k0-M*AxialWaveNumber_withFlow)^2;
            Res = k_y*tan(2*k_y*DuctHeight) - 1i*FluidProperties.Density*Omega*Y_ac *(1-M*AxialWaveNumber_withFlow/k0)^2;
        end     
        
        function [AxialWaveNumber_withFlow,AxialWaveNumber_AgainstFlow] = AxialWaveNumber(M,k,k_x,k_y)
            %Obtain the axial wavenumber in the duct, based on the
            %dispersion relation. 
            %The dispersion relation is obtained under the assumption of
            %uniform flow
            %There are two solutions, corresponding to
            %waves travelling with the flow and against the flow if M > 0. 
            %The positive x-direction is with the flow and thus if M = 0, 
            %the wave with the flow correspond to the waves propagating
            %in the positive axial direction.
            %Which solution that is taken depends on the sign of the
            %Mach-number            
            AxialWaveNumber_withFlow = k/(1-M^2)*(-M + sqrt(1 - (1-M^2)*(k_x^2/k^2+k_y^2/k^2)));
            AxialWaveNumber_AgainstFlow = k/(1-M^2)*(-M - sqrt(1 - (1-M^2)*(k_x^2/k^2+k_y^2/k^2)));

        end
    end
    
end

