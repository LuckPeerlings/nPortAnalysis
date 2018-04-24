classdef PronyMethod < handle
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
        %Vector of the measured pressures
        Frequency
        P
        MicSpacing %The microphone spacing
        NrModes
        AllAmplitudes
        WaveNumber
        AllWaveNumbers
        Epsilon
        Method
    end
    
    methods
        %Empty Constructor 
        %Constructor
        function obj = PronyMethod(Frequency,P,MicSpacing,NrModes,Epsilon,Method)
            if nargin > 0 
                obj.Frequency = Frequency; 
                obj.P = P;
                obj.MicSpacing = MicSpacing;
                obj.NrModes = NrModes;
                obj.Parameter = Parameter;
            elseif nargin == 1
                obj.TestFunctionNicolas
            else
                
                obj.TestClass;
            end
        end
        
        function obj = DetermineImpedance(obj)
            
        end
        
        function obj = DetermineWaveNumber(obj)
            %Calculating the wavenumber and the amplitude of the modes for
            %each frequency
            for ff = 1:size(obj.P,2)
                if strcmp(obj.Method, 'ESPRIT')
                    [Amplitude_temp,kappa_temp] = obj.ESPRIT(obj.P(:,ff).',obj.MicSpacing,obj.NrModes, obj.Epsilon); % 1e-10
                    obj.AllAmplitudes(:,ff) = [Amplitude_temp ; zeros(obj.NrModes - length(Amplitude_temp), 1)];
                    obj.AllWaveNumbers(:,ff) = [kappa_temp ; zeros(obj.NrModes - length(kappa_temp), 1)];
                elseif strcmp(obj.Method, 'MatrixPencil')
                    [Amplitude_temp,kappa_temp] = obj.MatrixPencil(obj.P(:,ff).',obj.MicSpacing,obj.NrModes, obj.Epsilon); % 1e-10
                    obj.AllAmplitudes(:,ff) = [Amplitude_temp ; zeros(obj.NrModes - length(Amplitude_temp), 1)];
                    obj.AllWaveNumbers(:,ff) = [kappa_temp ; zeros(obj.NrModes - length(kappa_temp), 1)];
                elseif strcmp(obj.Method, 'Basic')
                    [Amplitude_temp,kappa_temp] = obj.BasicPronyMethod(obj.P(:,ff).',obj.MicSpacing,obj.NrModes, obj.Epsilon);
                    obj.AllAmplitudes(:,ff) = [Amplitude_temp ; zeros(obj.NrModes - length(Amplitude_temp), 1)];
                    obj.AllWaveNumbers(:,ff) = [transpose(kappa_temp) ; zeros(obj.NrModes - length(kappa_temp), 1)];
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
            PronyMethod.ArrowPlotComplexDomain(AX,obj.kappa(nn,I))
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
        
        function TestClass(obj)
            % Set the properties of the class
            obj.NrModes = 12;
            obj.MicSpacing = 0.1;
            obj.Epsilon = 1-1e-10;
                                   
            %Function to determine whether the prony method is implemented
            %correctly
             %1-MicSpacing; %0.15
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            
            %The wavenumbers are be sorted using the size of the absolute
            %value.
            k = [2+1i,1+2i,5+4i,3,6i, 1.5+1.5i]*5;
            A = [10-8i,5+9i,-5-6i,-5+4i,2+2i,1-1i];
            
            %Determine pressure with the above exponentials at the
            %positions X.
            X = 0:obj.MicSpacing:2.4;
            Pressure = zeros(1,length(X));
            for ii = 1:length(k) %NrModes
                Pressure = Pressure + A(ii) * exp(1i * k(ii) * X);
            end
            
            obj.P = Pressure;
            
         
            [C, kappa] = PronyMethod.BasicPronyMethod(obj.P,obj.MicSpacing,obj.NrModes,obj.Epsilon);

            %Sort the modes by amplitude
            figure
            subplot(3,2,1)
            plot(X,real(Pressure));
            xlabel('Distance')
            ylabel('Real Pressure')

            subplot(3,2,2)
            plot(X,imag(Pressure));
            xlabel('Distance')
            ylabel('Imag Pressure')
            
            subplot(3,2,3)
            plot(real(kappa),'x')
            hold all
            plot(real(k),'o')
            ylabel('Real Wavenumber')
            xlabel('Mode nr')
            legend('Educed','Input')
            
            subplot(3,2,4)
            plot(imag(kappa),'x')
            hold all
            plot(imag(k),'o')
            ylabel('Imag Wavenumber')
            xlabel('Mode nr')   
            legend('Educed','Input')
                        
            subplot(3,2,5)
            plot(real(C),'x')
            hold all
            plot(real(A),'o')
            ylabel('Real Amplitude')
            xlabel('Mode nr')
            legend('Educed','Input')
            
            subplot(3,2,6)
            plot(imag(C),'x')
            hold all
            plot(imag(A),'o')
            ylabel('Imag Amplitude')
            xlabel('Mode nr')
            legend('Educed','Input')
        end
    end
    
    methods(Static)

        function [Amplitudes, kappa] = MatrixPencil(P, dX, L, epsilon)
            % L = number of modes that are calculated
            % M = number of modes we look at
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif L < 3
                error('The number of modes must be superior or equal to 3')
            end
            
            H = PronyMethod.Hankel(2*N-L, L+1, 0, P);
            
            [~, R, MPi] = qr(H);
            
            M = 0;
            energy = 0;
            Rabs = abs(diag(R));
            sum_Rabs = sum(Rabs);
            Rabs = Rabs / sum_Rabs;
            while M < length(Rabs) && energy < epsilon
                M = M+1;
                energy = energy + Rabs(M);
            end
            
            if M == 0
                M = 1;
            end
            
            S = R * transpose(MPi);
                        
            T_s0 = PronyMethod.Matrix_T(S, M, L, 0);
            T_s1 = PronyMethod.Matrix_T(S, M, L, 1);
            
            D = diag(diag(R(1:M, 1:M)));
            
            T_s0 = pinv(D) * T_s0;
            T_s1 = pinv(D) * T_s1;
            
            F = pinv(transpose(T_s0)) * transpose(T_s1);
            
            Z = eig(F);
            
            if length(Z) < M
                disp('There are less eigenvalues than modes')
            elseif length(Z) > M
                disp('There are more eigenvalues than modes')
            end
            
            f_i = log(Z);
            
            V = PronyMethod.Vandermonde(2*N, Z);
            
            C = V \ transpose(P);
            
            Amplitudes = C;
            [~,I] = sort(abs(C),'descend');
            Amplitudes = Amplitudes(I);
            
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            kappa = kappa(I);
            kappa = kappa.';
            
        end
        
        function [Amplitudes, kappa] = BasicPronyMethod(P, dX, L, epsilon)
            % L = number of modes that are calculated
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            end
            
            H = PronyMethod.Hankel(2*N-L, L, 0, P);
            h = P(L+1:2*N).';
            
%             H_2 = toeplitz(P(L:(2*N-1)),fliplr(P(1:L)));
%             h_2 = P(L+1:end).';
            
            Z_0 = roots([1 ; flip(-pinv(H)*h)]);
            
            V_0 = PronyMethod.Vandermonde(2*N, Z_0);
            
            C_0 = V_0\transpose(P);
            
            index = [];
            for ii = 1:length(C_0)
                if abs(C_0(ii)) > 1-epsilon
                    index = [index, ii];
                end
            end
            
            Z_1 = Z_0(index);
            
            V_1 = PronyMethod.Vandermonde(2*N, Z_1);
            
            C_1 = V_1\transpose(P);
            
            f_i = log(Z_1);
            
            Amplitudes = C_1;
            [~,I] = sort(abs(C_1),'descend');
            Amplitudes = Amplitudes(I);
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            kappa = kappa(I);
            
        end
        
        function [Amplitudes, kappa] = ESPRIT(P, dX, L, epsilon)
            % L = number of modes that are calculated
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            end
            
            H = PronyMethod.Hankel(2*N-L, L+1, 0, P);
            
            [~, D, W] = svd(H);
            
            W = ctranspose(W);
            
            sigma = diag(D);
            
            sigma = abs(sigma);
%             sigma = sigma.^2;
            sum_sigma = sum(sigma);
            sigma = sigma/sum_sigma;
            energy = 0;
            
            M = 0;
            while M < length(sigma) && energy < epsilon
                M = M + 1;
                energy = energy + sigma(M);
            end
            
            W_rect = W(1:M, 1:L+1);
            
            W_0 = W_rect(1:M, 1:L);
            W_1 = W_rect(1:M, 2:L+1);
            
            F = pinv(transpose(W_0)) * transpose(W_1);
            
            Z = eig(F);
            
            V = PronyMethod.Vandermonde(2*N, Z);
            
            C = V\transpose(P);
            
            f_i = log(Z);
            
            kappa = -1i * f_i / dX;
            
            Amplitudes = C;
            [~,I] = sort(abs(C),'descend');
            Amplitudes = Amplitudes(I);
            
            kappa = kappa(I);
            
        end
        
        function [AX] = ArrowPlotComplexDomain(AX,Data)
            
            for ii = 1:length(Data)-1
                if ii == 1
                quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'o','filled','MarkerFaceColor','r','MaxHeadSize',0.5);
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
    end
    
end
