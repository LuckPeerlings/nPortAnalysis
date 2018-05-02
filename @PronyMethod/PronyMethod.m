classdef PronyMethod < matlab.mixin.SetGet
    %This class is used to determine the impedance using the so called
    %Prony method.
    %The details of this method are taken from the paper
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %TODO List
    % - Normalization or any other method to be able to define epsilon in a
    %   logical matter *check* matrix norm ? sing value norm. ? *check*
    % - Fix the MyLittleProny *check*
    % - Fix the code so that it works for inputs of P which are matrices
    % *check*
    % - Include code to be able to get the data on regular intervals
    % - Include a comparison of the different methods with the inclusion of
    %   noise on the input data.
    % - Comment the code with references to the paper where necessary
    
    properties
        Frequency      % Column vector of the frequencies used during the measurement, in increasing order 
        P              % Matrix of the measured complex pressures, of size (number of frequencies used) * (number of measurement points)
        EqP            % Pressures on equispaced grid
        
        NrModes        % Integer, maximum number of modes we want to calculate        
        WaveNumber     % Row vector of the wavenumber corresponding to the first mode for each frequency used       
        Epsilon        % Real number between 0 and 1, see description of its use in BasicPronyMethod, MatrixPencil and ESPRIT
        
        Method = 'ESPRIT';         % Character string, name of the method used to determine the amplitudes and wavenumbers of each mode of the complex signal (choice between 'BasicPronyMethod', 'MatrixPencil' and 'ESPRIT')
        AllWaveNumbers % Matrix of the wavenumbers of each mode and for each frequency used, of size (NrModes) * (number of frequencies used)
        AllAmplitudes  % Matrix of the amplitudes of each mode and for each frequency used, of size (NrModes) * (number of frequencies used) 
        
        %Private properties
        MicPositions   % Microphone positions used in the measurement
        MicEqPositions % Equidistant microphone positions used to evaluate the prony method
        MicSpacing     % Real number, distance between the measurement points
    end
    
    methods
        %Empty Constructor 
        %Constructor
        function obj = PronyMethod(Frequency,P,MicPositions,MicEqPositions,NrModes,Epsilon)
            if nargin > 1 
                obj.Frequency = Frequency; 
                obj.P = P;
                obj.MicPositions = MicPositions;
                obj.MicEqPositions = MicEqPositions;
                obj.NrModes = NrModes;
                obj.Epsilon = Epsilon; %1-1e-10; % This value of Epsilon is found to function best.
                if length(uniquetol(diff(MicEqPositions)))>1
                    error('The given equidistant microphone positions are not equidistant')
                end
                obj.MicSpacing = obj.MicEqPositions(2)-MicEqPositions(1);
            elseif nargin == 1
                obj.TestFunctionNicolas
            else                
                %obj.TestClass;
            end
        end
        function obj = EquispaceData(obj)
            for ii = 1:size(obj.P,2)
                [EquiPressure(:,ii)] = PronyMethod.equispacing(obj.MicPositions,obj.MicEqPositions,obj.P(:,ii),length(obj.P(:,ii)),12.3);
            end
            obj.EqP = EquiPressure;
        end
        function obj = DetermineWaveNumber(obj,Method)
            %Calculating the wavenumber and the amplitude of the modes for
            %each frequency
            %The default method is the ESPRIT method.
            if nargin == 1    
                Method = obj.Method;                
            end
            
            for ff = 1:size(obj.EqP,2)
                if strcmp(Method, 'ESPRIT')
                    [obj.AllAmplitudes(:,ff),obj.AllWaveNumbers(:,ff)] = obj.ESPRIT(obj.EqP(:,ff).',obj.MicSpacing,obj.NrModes, obj.Epsilon); % 1e-10
                elseif strcmp(Method, 'MatrixPencil')
                    [obj.AllAmplitudes(:,ff),obj.AllWaveNumbers(:,ff)] = obj.MatrixPencil(obj.EqP(:,ff).',obj.MicSpacing,obj.NrModes, obj.Epsilon); % 1e-10
                elseif strcmp(Method, 'Basic')
                    [obj.AllAmplitudes(:,ff),obj.AllWaveNumbers(:,ff)] = obj.BasicPronyMethod(obj.EqP(:,ff).',obj.MicSpacing,obj.NrModes, obj.Epsilon);
                end
            end
            obj.WaveNumber = obj.AllWaveNumbers(1,:);
        end
        
        function obj = PlotWaveNumberComplexDomain(obj,ModeNr)
            nn = ModeNr;
            
            [~,I] = sort(obj.Frequency);
            figure;
            AX  = axes();
            hold on;
            PronyMethod.ArrowPlotComplexDomain(AX,obj.AllWaveNumbers(nn,I))
            axis square;
            
            ylabel('Imag. part of wavenumber')  
            xlabel('Real part of wavenumber')
           
        end
        
        function obj = PlotAmplitude(obj)
            %Function to plot the obtained real and imaginary part of each
            %of the modes and the relative amplitude of the modes compared to the largest mode present. 
            
            figure
            for nn = 1:obj.NrModes
                subplot(obj.NrModes+1,1,nn)
                plot(obj.Frequency,real(obj.kappa(nn,:)),'x');
                hold all
                plot(obj.Frequency,imag(obj.kappa(nn,:)),'o');
            end
            legend('Real part of wavenumber','Imag. part of wavenumber')  
            
           
            AX1 = subplot(obj.NrModes+1,1,obj.NrModes+1);            
            for ff = 1:size(obj.Frequency,2)
                plot(obj.Frequency(ff),abs(obj.Amplitude(:,ff))/norm(obj.Amplitude(:,ff)),'x')
                hold on         
                AX1.ColorOrderIndex = 1;
                ylabel('Rel. Amplitude')
                xlabel('Frequency')
            end
        end
        
        function TestFunctionNicolas(obj)
        %function to determine whether the equispacing function works well
        
        %different signals to be tested 
        obj.MicEqPositions = 0:0.05:5*0.05;
        obj.MicSpacing = obj.MicEqPositions(2)-obj.MicEqPositions(1);
            
        obj.MicPositions = obj.MicEqPositions;
        obj.MicPositions(2:end-1) = obj.MicEqPositions(2:end-1) + 0.01*obj.MicSpacing/2*rand(1,length(obj.MicPositions)-2); 
        
        X_test=transpose(obj.MicPositions);
        t= transpose(0:0.01:1);
        equiX= transpose(obj.MicEqPositions);
        equigrid=t;
        
        func = @(x) exp(1i*(-11.510242292516846 + 6.292327191264097i)*x);
        Testsignal=func(X_test);
        Testsignal2=func(t);
        %function call: 
        % 1. Parameter: KxN Vector of space locations
        % 2. Parameter: KxN Vector of testfunction
        % 3. Parameter: Value of L: K>=L>=N with N= nonharmonic bandwidth
        % 4. Parameter: Vaue of b:  Variance of the gaussian
        
        [equiPressure] =PronyMethod.equispacing(X_test,equiX,Testsignal,6,12.3);
        [Testdata] =PronyMethod.equispacing(t,equigrid,Testsignal2,99,12.3);
        
%        MSE of the approximated Data and the Testdata
         error=equiPressure-func(equiX);

        figure
        subplot(3,2,1)
        plot(equigrid,real(Testdata),'x')
        ylabel('Reconstructed Amplitude')
        xlabel('Equispaced Steps')
        hold on
        plot(t,real(Testsignal2))
        plot(equiX,real(equiPressure),'o')
        title('Real part')
      
        
        %imaginary Part
        
        subplot(3,2,2)
        plot(equigrid,imag(Testdata),'x')
        ylabel('Reconstructed')
        xlabel('Equispaced Steps')
        hold on
        plot(t,imag(Testsignal2))
        plot(equiX,imag(equiPressure),'o')
        title('Imaginary part')
        legend('generated Data','Testsignal','Measurement points')
        
        
        % Plot errrors of real Part  
        subplot(3,2,3)
        plot(equiX,real(error))
        ylabel('Error of reconstruction')
        xlabel('Equispaced Steps')
        hold on
        plot(t,zeros(length(t)));
        
        % Plot errrors of imag Part  
        subplot(3,2,4)
        plot(equiX,imag(error))
        ylabel('Error of reconstruction')
        xlabel('Equispaced Steps')        
        hold on
        plot(t,zeros(length(t)));
        
        % plot of the input with non equispaced nodes 
        %real
        
        subplot(3,2,5)
        plot(X_test,real(Testsignal),'o')
        ylabel('Testfunction')
        xlabel('Non-Equispaced Steps')
        hold on
        plot(t,real(Testsignal2))
        
        %imaginary
        subplot(3,2,6)
        plot(X_test,imag(Testsignal),'o')
        ylabel('Testfunction')
        xlabel('Non-Equispaced Steps')
        hold on
        plot(t,imag(Testsignal2))
        
        %best parameter estimate for b
        kx=[4.48291066436601 - 0.0232247424295838i;5.41782308166276 - 0.0252505740248546i;6.34711835706620 - 0.0271515381332884i;7.27292078959503 - 0.0289422972352744i;8.19639320338978 - 0.0306373246438771i;9.11822807268444 - 0.0322490790633878i;10.0388646289507 - 0.0337879280000863i;10.9585955258497 - 0.0352624496801973i;11.8776236209035 - 0.0366797698925568i;12.7960941782502 - 0.0380458519566265i;13.7141140953456 - 0.0393657289147092i;14.6317638826141 - 0.0406436852143972i;15.5491054014463 - 0.0418833981827764i;16.4661870190230 - 0.0430880485163045i;17.3830471353964 - 0.0442604071291219i;18.2997166540144 - 0.0454029039623454i;19.2162207483640 - 0.0465176829647735i;20.1325801487590 - 0.0476066463992921i;21.0488120952079 - 0.0486714908444128i;21.9649310535897 - 0.0497137366808734i;22.8809492612396 - 0.0507347524250974i;23.7968771477129 - 0.0517357749536686i;24.7127236629481 - 0.0527179264258585i;25.6284965358633 - 0.0536822285329635i;26.5442024800763 - 0.0546296145681435i;27.4598473590062 - 0.0555609397073330i;28.3754363194589 - 0.0564769898124638i;29.2909739005355 - 0.0573784890067363i;30.2064641230508 - 0.0582661062236577i;31.1219105634366 - 0.0591404608938080i;32.0373164152010 - 0.0600021279033990i;32.9526845403354 - 0.0608516419348758i;33.8680175125510 - 0.0616895012807165i;34.7833176538282 - 0.0625161712061832i;35.6985870654663 - 0.0633320869242929i;36.6138276545808 - 0.0641376562360900i];
        b=linspace(1,100,100);
        nn=length(kx);
        predicted=zeros(length(X_test),length(b),nn);
        real1=zeros(length(X_test),length(b),nn);
        fonc = @(x,y) exp(1i*(y)*x);
        
        
        for bb=1:length(b)
            for kk=1:nn
            [predicted(:,bb,kk)]=PronyMethod.equispacing(X_test,equiX,fonc(X_test,kx(kk)),length(X_test),b(bb));
            real1(:,bb,kk)=fonc(X_test,kx(kk));
            end
        end
        argument=predicted-real1;
        handover=sqrt(dot(argument,argument));
        mse=(1/nn)*dot(sqrt(handover),sqrt(handover),3);
        
        figure
        plot(b,mse)
        end
        
        function obj = TestClass(obj)
            % Set the properties of the class
            
            obj.Epsilon = 15e-2; 1e-10; %1-1e-10;  
            
            
            obj.MicEqPositions = 0:0.1:1;
            obj.MicSpacing = obj.MicEqPositions(2)-obj.MicEqPositions(1);
            
            obj.MicPositions = obj.MicEqPositions;
            obj.MicPositions(2:end-1) = obj.MicEqPositions(2:end-1) + 1*obj.MicSpacing/2*rand(1,length(obj.MicPositions)-2);
            
            obj.NrModes = floor(length(obj.MicPositions)/2)-3;
            
            %Function to determine whether the prony method is implemented
            %correctly
            %1-MicSpacing; %0.15
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            
            %The wavenumbers are be sorted using the size of the absolute
            %value.
            k = [3-0.3i,3-0.2i,15-0.12i,9, 4.5-0.45i];
            A = [10-10i,5+10i,-5-6i,-5+4i,2+2i,1-1i];
%             
              
            k = k(1:1);
            A = A(1:1);
            
            %Determine pressure with the above exponentials at the
            %positions X.            
            Pressure = zeros(1,length(obj.MicPositions));
            for ii = 1:length(k) %NrModes
                Pressure = Pressure + A(ii) * exp(1i * k(ii) * obj.MicPositions);
            end
            SNR = 0.05;
            L = length(Pressure);
            E = mean((abs(Pressure)).^2);
            noise = sqrt(SNR * E / 2) .* (randn(1,L) + 1i.*randn(1,L));
            obj.P = Pressure + noise;
            size(obj.P)
            isrow(obj.P)
            [equiPressure] = PronyMethod.equispacing(obj.MicPositions.',obj.MicEqPositions.',obj.P.',length(obj.P),12.3);
            %Approximate the function values on a randomized non-equispaced grid
           
            
                        

            [C, kappa] = PronyMethod.BasicPronyMethod(equiPressure.',obj.MicSpacing,obj.NrModes,obj.Epsilon);
            [C, kappa] = PronyMethod.ESPRIT(equiPressure.',obj.MicSpacing,obj.NrModes,obj.Epsilon);

            %Sort the modes by amplitude
            figure
            subplot(3,2,1)
            plot(obj.MicPositions,real(Pressure),'o');
            hold all
            plot(obj.MicEqPositions,real(equiPressure),'x');
            xlabel('Distance')
            ylabel('Real Pressure')

            subplot(3,2,2)
            plot(obj.MicPositions,imag(Pressure),'o');
            hold all
            plot(obj.MicEqPositions,imag(equiPressure),'x');
            xlabel('Distance')
            ylabel('Imag Pressure')
            
            subplot(3,2,3)
            plot(real(kappa),'x')
            hold all
            plot(real(k),'o')
            ylabel('Real Wavenumber')
            xlabel('Mode nr')
%             legend('Educed','Input')
            
            subplot(3,2,4)
            plot(imag(kappa),'x')
            hold all
            plot(imag(k),'o')
            ylabel('Imag Wavenumber')
            xlabel('Mode nr')   
%             legend('Educed','Input')
                        
            subplot(3,2,5)
            plot(real(C),'x')
            hold all
            plot(real(A),'o')
            ylabel('Real Amplitude')
            xlabel('Mode nr')
%             legend('Educed','Input')
            
            subplot(3,2,6)
            plot(imag(C),'x')
            hold all
            plot(imag(A),'o')
            ylabel('Imag Amplitude')
            xlabel('Mode nr')
%             legend('Educed','Input')
        end
    end
    
    methods(Static)
        
        % The three following fonctions are all different methods to
        % calculate the modeshapes of a measured signal. All of these
        % methods are based on the assumption that the signal can be
        % expressed as a sum of damped exponentials. Thus, the goal of
        % these methods is to calculate the wavenumbers and the
        % corresponding amplitudes of these exponentials for a certain
        % number of modes, number which is determined inside the
        % functions. All of these functions are based on the so-called
        % "Prony Method", the first function corresponding to the original
        % method, and the all of the algorithms on which their implementations are
        % based come from [D. Potts, M. Tasche/Linear Algebra and its
        % Applications 439 (2013) 1024-1039].
        
        function [Amplitudes, Wavenumbers] = BasicPronyMethod(P, dX, L, epsilon)
            % This method is one of the three Prony based methods designed to compute the wavenumbers and amplitudes of each mode of the Prony decomposition of a measured signal.
            %
            % Inputs :
            % - P = measured signal (Pressure)
            % - dX = distance between the measurement points (MicSpacing)
            % - L = number of modes that are calculated (upper bound for M)
            % - epsilon = criterion used to discriminate low-influence modes
            %
            % Outputs :
            % - Amplitudes : amplitudes corresponding to each guessed mode (column vector)
            % - kappa : wavenumbers corresponding to each guessed mode (column vector)
            
            % The principle of this method is explained from page 1027 to
            % 1029 of [D. Potts, M. Tasche] (see reference under "methods 
            % (Static)", above). The algorithm is displayed at page 1029. 
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif epsilon < 0 || epsilon > 1
                error('Epsilon must be a number between 0 and 1')
            end
            
            % First, the exponentials have to be calculated. These
            % exponentials can be deduced from the roots of q, the "Prony 
            % polynomial", whose coefficients vector is the solution of the
            % Yule-Walker system H*q = h.
            
            % Thus, in the first place, the Yule-Walker system is built :
            
            H = PronyMethod.Hankel(2*N-L, L, 0, P);
            h = P(L+1:2*N).';
            
            % Then, the roots of the polynomial built with its solution are
            % calculated :
            
            Z_0 = roots([1 ; flip(-pinv(H)*h)]);
            
            % The roots of the polynomial are assumed to be the 
            % exponentials. Now, the corresponding amplitudes have to be 
            % calculated.
            
            % As the measurement points are equally spaced, one can assume
            % that the coefficients of the exponentials are the solution of
            % the system V*C = P, where P is the vector of the measured
            % data and V a Vandermonde matrix based on the vector of the
            % exponentials (as powering the exponentials would be
            % equivalent to incrementing the measurement point).
            
            % First the Vandermonde is built via the "Vandermonde"
            % function (implemented after the three methods) :
            
            V_0 = PronyMethod.Vandermonde(2*N, Z_0);
            
            % Then the amplitudes are calculated :
            
            C_0 = V_0\transpose(P);
            
            % Now, we have a set of modes, whose wavenumbers can be
            % easily deduced as the distance between the measurement points
            % is known, and their corresponding amplitudes. The next step 
            % is to get rid of the modes whose amplitude are low compared
            % to the other ones.
            
            % In order to discriminate the lowest modes, the vector is
            % sorted (in decreasing order), and a variable "energy",
            % corresponding to a proportion of the value of the sum of all
            % of the amplitudes. The highest amplitudes (in terms of 
            % absolute value) are chosen, until the sum of the chosen 
            % amplitudes reaches a certain amount of energy (epsilon). 
            % This very operation is used in the three methods.
            
            C_0_abs = abs(C_0);
            sum_C_0 = sum(C_0_abs);
            C_0_abs = C_0_abs / sum_C_0;
            C_0_abs = sort(C_0_abs, 'descend');
            energy = 0;
            index = 0;
            while index < length(C_0_abs) && energy < epsilon
                index = index+1;
                energy = energy + C_0_abs(index);
            end
            
            % Safety case :
            
            if index == 0
                index = 1;
            end
            
            % Now that we know what modes have the most influence in the 
            % signal (in terms of amplitude), we only keep their
            % corresponding exponentials :
            
            Z_1 = Z_0(1:index);
            
            % A new Vandermonde matrix is then built with the new
            % exponential vector :
            
            V_1 = PronyMethod.Vandermonde(2*N, Z_1);
            
            % The new amplitudes are calculated :
            
            C_1 = V_1\transpose(P);
            
            % The wavenumbers (kappa) can now be deduced from the
            % exponentials :
            
            f_i = log(Z_1);
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            % At last, the modes are sorted by decreasing amplitudes :
            
            
            Amplitudes = zeros(1,L);
            [~,I] = sort(abs(C_1),'descend');
            Amplitudes(1:length(I)) = C_1(I);
            
            Wavenumbers = zeros(1,L);
            Wavenumbers(1:length(I)) = kappa(I);
        end
        
        function [Amplitudes, Wavenumbers] = MatrixPencil(P, dX, L, epsilon)
            % This method is one of the three Prony based methods designed to compute the wavenumbers and amplitudes of each mode of the Prony decomposition of a measured signal.
            %
            % Inputs :
            % - P = measured signal (Pressure)
            % - dX = distance between the measurement points (MicSpacing)
            % - L = number of modes that are calculated (upper bound for M)
            % - epsilon = criterion used to discriminate low-influence modes
            %
            % Outputs :
            % - Amplitudes : amplitudes corresponding to each guessed mode (column vector)
            % - kappa : wavenumbers corresponding to each guessed mode (column vector)
            
            % The principle of this method is explained oat pages 1030 and
            % 1031 of [D. Potts, M. Tasche] (see reference under "methods 
            % (Static)", above). The algorithm is displayed at page 1032.
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif L < 3
                error('The number of modes must be superior or equal to 3')
            elseif epsilon < 0 || epsilon > 1
                error('Epsilon must be a number between 0 and 1')
            end
            
            % In this method, we try to solve the generalized eigenvalues
            % problem z*H(0) - H(1), where H(0) and H(1) are Hankel 
            % matrices of measured data (see the "Hankel" function after 
            % the three methods). As H is, in general, rank deficient, we
            % use a qr-decomposition to evaluate its order and simplify the
            % problem.
            
            % First, the H matrix is built :
            
            H = PronyMethod.Hankel(2*N-L, L+1, 0, P);
            
            % Then, a qr-decomposition is made on the H matrix. The
            % matrices R and MPi (transposition matrix) are the only ones
            % interesting.
            
            [~, R, MPi] = qr(H);
            
            % In order to discriminate the lowest modes, a variable
            % "energy", corresponding to a proportion of the value of the 
            % sum of all of the diagonal elements of R, is created. The 
            % highest elements (in terms of absolute value) are chosen, 
            % until the sum of the chosen elements reaches a certain 
            % amount of energy (epsilon). This very operation is used in
            % the three methods.
            
            M = 0;
            energy = 0;
            Rabs = abs(diag(R));
            sum_Rabs = sum(Rabs);
            Rabs = Rabs / sum_Rabs;
            while M < length(Rabs) && energy < epsilon
                M = M+1;
                energy = energy + Rabs(M);
            end
            
            % Safety case :
            
            if M == 0
                M = 1;
            end
            
            % Now that we have found how many modes we want to keep, we can
            % re-build the problem in order to find its solutions, with
            % H = Q * S :
            
            S = R * transpose(MPi);
            
            % As Q is unitary, the solutions z of our problem are the same
            % as for the problem z*S(0) - S(1), which are the same as  for
            % the problem z*T(0) - T(1), T being equal to S without its
            % null values (from index M+1 to 2*N-L-M), we re-build the
            % problem with matrices T(0) and T(1) :
                        
            T_s0 = PronyMethod.Matrix_T(S, M, L, 0);
            T_s1 = PronyMethod.Matrix_T(S, M, L, 1);
            
            % As the first M diagonal elements of R can be used as a
            % preconditioning matrix, we build the corresponding matrix :
            
            D = diag(diag(R(1:M, 1:M)));
            
            % Then the matrices T are conditioned :
            
            T_s0 = pinv(D) * T_s0;
            T_s1 = pinv(D) * T_s1;
            
            % The generalized eigenvalues problem can then be transformed
            % into a simple problem :
            
            F = pinv(transpose(T_s0)) * transpose(T_s1);
            
            % The exponentials can now be computed as the eigenvalues of
            % matrix F :
            
            Z = eig(F);
            
            if length(Z) < M
                disp('There are less eigenvalues than modes')
            elseif length(Z) > M
                disp('There are more eigenvalues than modes')
            end
            
            % As the measurement points are equally spaced, one can assume
            % that the coefficients of the exponentials are the solution of
            % the system V*C = P, where P is the vector of the measured
            % data and V a Vandermonde matrix based on the vector of the
            % exponentials (as powering the exponentials would be
            % equivalent to incrementing the measurement point).
            
            % First the Vandermonde is built via the "Vandermonde"
            % function (implemented after the three methods) :
            
            V = PronyMethod.Vandermonde(2*N, Z);
            
            % The amplitudes are now calculated :
            
            C = V \ transpose(P);
            
            % The wavenumbers (kappa) can now be deduced from the
            % exponentials :
            
            f_i = log(Z);
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            % At last, the modes are sorted by decreasing amplitudes :
            
            Amplitudes = zeros(1,L);
            [~,I] = sort(abs(C),'descend');
            Amplitudes(1:length(I)) = C(I);
            
            Wavenumbers = zeros(1,L);
            Wavenumbers(1:length(I)) = kappa(I);
        end
        
        function [Amplitudes, Wavenumbers] = ESPRIT(P, dX, L, epsilon)
            % This method is one of the three Prony based methods designed to compute the wavenumbers and amplitudes of each mode of the Prony decomposition of a measured signal.
            %
            % Inputs :
            % - P = measured signal (Pressure)
            % - dX = distance between the measurement points (MicSpacing)
            % - L = number of modes that are calculated (upper bound for M)
            % - epsilon = criterion used to discriminate low-influence modes
            %
            % Outputs :
            % - Amplitudes : amplitudes corresponding to each guessed mode (column vector)
            % - kappa : wavenumbers corresponding to each guessed mode (column vector)
            
            % The principle of this method is explained on pages 1032 and
            % 1033 of [D. Potts, M. Tasche] (see reference under "methods 
            % (Static)", above). The algorithm is displayed at page 1033.
            if ~isrow(P)
                error('Pressure vector has to be a row vector')
            end
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif epsilon < 0 || epsilon > 1
                error('Epsilon must be a number between 0 and 1')
            end
            
            % The principle of the ESPRIT method is quite the same as for
            % the Matrix Pencil method, but the qr-decomposition is
            % replaced by a singular values decomposition (SVD). Thus, the
            % purpose is to calculate the exponentials of the Prony
            % decomposition as the solution of the generalized eigenvalues
            % problem z*H(0) - H(1) where H is a Hankel matrix of measured
            % data.
            
            % In the first place, the matrix H is built :
            
            H = PronyMethod.Hankel(2*N-L, L+1, 0, P);
            
            % Then, a SVD is operated on the matrix H :
            
            [~, D, W] = svd(H);
            
            W = ctranspose(W);
            
            % Truncate the SVD such that sigma(M+1) < epsilon*sigma(1),
            % Rebuild the matrix with M modes.
            
            sigma = diag(D);
            M = 1;
            while sigma(M) >= epsilon*sigma(1) && M < length(sigma)
                M = M+1;

            end
            M = M-1;
            
                                    
            % Now, the problem can be written with only the selected modes.
            % For this, the matrix W has to be re-written :
            
            W_rect = W(1:M, 1:L+1);
            
            % In order to correspond to the matrices H(0) and H(1), we
            % define analog matrices W(0) and W(1) so that H(s) = U*D*W(s).
            
            W_0 = W_rect(1:M, 1:L);
            W_1 = W_rect(1:M, 2:L+1);
            
            % As U is unitary, the problem can be written z*D*W(0)-D*W(1),
            % and as one can multiply the problem with D^(-1), the problem
            % can be written z*transpose(W(0))-transpose(W(1)).
            
            % The problem can then be transformed into a simple eigenvalue
            % problem of the matrix F defined by :
            
            F = pinv(transpose(W_0)) * transpose(W_1);
            
            % The exponentials can be computed as the eigenvalues of F :
            
            Z = eig(F);
            
            % The wavenumbers can the be deduced from the exponentials :
            
            f_i = log(Z);
            kappa = -1i * f_i / dX;
            
            % As the measurement points are equally spaced, one can assume
            % that the coefficients of the exponentials are the solution of
            % the system V*C = P, where P is the vector of the measured
            % data and V a Vandermonde matrix based on the vector of the
            % exponentials (as powering the exponentials would be
            % equivalent to incrementing the measurement point).
            
            % First the Vandermonde is built via the "Vandermonde"
            % function (implemented after the three methods) :
            
            V = PronyMethod.Vandermonde(2*N, Z);
            
            % The amplitude can then be computed :
            
            C = V\transpose(P);
            
            % At last, the modes are sorted by decreasing amplitudes :
            
            Amplitudes = zeros(1,L);
            [~,I] = sort(abs(C),'descend');
            Amplitudes(1:length(I)) = C(I); 
            
            Wavenumbers = zeros(1,L);
            Wavenumbers(1:length(I)) = kappa(I); 
        end
        
        function [AX] = ArrowPlotComplexDomain(AX,Data)
            
            for ii = 1:length(Data)-1
                if ii == 1
                quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'ko','filled','MarkerFaceColor','r','MaxHeadSize',0.5);
                else
                                    quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'ko','filled','MaxHeadSize',0.5);
                end                
            end
            plot(real(Data(ii+1)),imag(Data(ii+1)),'rx');
                       
        end
        
        function [H] = Hankel(M, N, s, P)
            H = zeros(M, N);
            for ii = 1:M
                for jj = 1:N
                    H(ii, jj) = P((ii - 1) + (jj - 1) + s + 1);
                end
            end
        end
        
        function [V] = Vandermonde(M, Z)
            for ii = 1:M 
                   V(ii, :) = Z(:) .^ (ii - 1);
            end
        end
        
        function [T] = Matrix_T(S, M, L, s)
           T = S(1:M, (1+s):(L+s));
        end
        
        function [S_1] = Matrix_S(S_0, L, s)
            S_1 = S_0(1:length(S_0(:,1)), (1+s):(L+s));
        end
        
       function [equiPressure] = equispacing(X, equiX, P,n,b)
         %This function approximates the complex pressure function with
         %nonequispaced steps and transforms it to a equispaced stepsized
         
         % Input 
         % 1. Parameter: Kx1 Vector of space locations, K= number nodes 
         % 2. Parameter: Kx1 Vector of equispace locations
         % 3. Parameter: Kx1 Vector of testfunction
         % 4. Parameter: Value of L: K>=L>=N with N= nonharmonic bandwidth 
         % 5. Parameter: b= Variance of the gaussian
        
         % Output
         % 1. Parameter: Kx1 vector equiPressure= approximated Pressure with homogenious
         %               step size 
         
           
         % 1) Initialization of parameters for the algorithm
         N=length(X);
                 
         % error exceptions to get the right dimensions 
         if  isrow(X)
             error('X has to be column vector');
         end
         
         if  isrow(equiX)
             error('equiX has to be column vector');
         end
         if  isrow(P)
             error('P has to be column vector');
         end
         
             
         %Preparing the window function and the stepsize     
         argumentX=zeros(N,n);
         augmentedX=zeros(N,n);
         
         l=-n/2:(n/2)-1;
         
         %argumentX is going to be the new argument used for the
         %windowfunction with eiquispaced input
         for ll=1:n
            argumentX(:,ll)=((X-max(X)/2)/max(X))-((1/n)*l(ll)); 
         end
         
         %augmented X is going to be the new argument used for the
         %windowfunction with eiquispaced input
         for ll=1:n
            augmentedX(:,ll)=((equiX-max(X)/2)/max(X))-((1/n)*l(ll)); 
         end         
         
         %truncated window function 
         truncatedphi= PronyMethod.twindow(argumentX,b,n);
         truncatedphi2= PronyMethod.twindow(augmentedX,b,n);

         % Computation of the support coefficients h
         h=truncatedphi\P; 

        % Computation of the pressure with equispaced nodes
         equiPressure = truncatedphi2*h;

        end
        
        function truncatedphi= twindow(argument,b,n)
            % this functions sets up the truncated window function for a
            % gaussian bell window function 
            % inputs: argument:  NxL Matrix of different locations in space per node
            %         b: Variance of the gaussian
            %         n: number of parameters h
           
        
            % Initialization of parameters
            truncatedphi=zeros(size(argument,1),n);  

            %truncated phi is a matrix Nxn  containing the different values
            %of the gaussian bell depending on the parameters b,n
            for ii=1:size(argument,1)
                for jj=1:n
                    truncatedphi(ii,jj)=exp(-((n*argument(ii,jj))^2)/b);
                end
            end
         end

    end
    
end
