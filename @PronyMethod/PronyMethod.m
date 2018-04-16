classdef PronyMethod < handle
    %This class is used to determine the impedance using the so called
    %Prony method.
    %The details of this method are taken from the paper
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Vector of the measured pressures
        Frequency
        P
        dX %The microphone spacing
        NrModes
        Amplitude
        kappa
        Parameter
    end
    
    methods
        %Constructor
        function obj = PronyMethod(Frequency,P,dX,NrModes,Parameter,noise)
            obj.Frequency = Frequency; 
            obj.P = P;
            obj.dX = dX;
            obj.NrModes = NrModes;
            obj.Parameter = Parameter;
        end
        
        function obj = DetermineImpedance(obj)
            
        end
        
        function obj = DetermineWaveNumber(obj)
            %Calculating the wavenumber and the amplitude of the modes for
            %each frequency
            for ff = 1:size(obj.P,2)
                if strcmp(obj.Parameter, 'Normal') == 1
                    [obj.Amplitude(:,ff),obj.kappa(:,ff)] = obj.PronyCoefficients(obj.P(:,ff).',obj.dX,obj.NrModes);
                elseif strcmp(obj.Parameter, 'Improved') == 1
                    [obj.Amplitude(:,ff),obj.kappa(:,ff)] = obj.PronyCoefficientsImproved(obj.P(:,ff).',obj.dX,obj.NrModes);
                end
            end
        end
        
        function obj = PlotWaveNumberComplexDomain(obj,ModeNr)
            nn = ModeNr;
            
            [c,I] = sort(obj.Frequency);
            I;
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
            %Function to determine whether the prony method is implemented
            %correctly
            NrModes = 20;
            MicSpacing = 1;
            X = 0:MicSpacing:39; %1-MicSpacing; %0.15
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            
            %The wavenumbers are be sorted using the size of the absolute
            %value.
            f = [7, 21, 200, 201, 53, 1000] * (1i / 1000); %[2+1i,1+2i,5+4i,3,6i, 1.5+1.5i]*10;
            A = [6,5,4,3,2,1]; %[15+9i,10-8i,2+2i,1-1i,-5+4i,-6-6i]; %,2+2i,1-1i
%             z = [0.9856-0.1628i, 0.9856+0.1628i, 0.8976-0.4305i, 0.8976+0.4305i, 0.8127-0.5690, 0.8127+0.5690];
            k = -1i * f / MicSpacing;
            
            noise = PronyMethod.Noise(length(X),(0.5+0.3i));
            
            %Determine pressure with the above exponentials
            Pressure = zeros(1,length(X));
            for ii = 1:length(k) %NrModes
                Pressure = Pressure + A(ii) * exp(1i * k(ii) * X);
            end
%             Pressure = Pressure + noise;
            
            if strcmp(obj.Parameter, 'Normal')
%                 [C, kappa] = obj.PronyCoefficients(Pressure,MicSpacing,NrModes);
%                 [C, kappa] = obj.MatrixPencil(Pressure, MicSpacing, NrModes);
%                 [C, kappa] = obj.MyLittleProny(Pressure, MicSpacing, NrModes, 1e-10);
                [C, kappa] = obj.ESPRIT(Pressure, MicSpacing, NrModes, 1e-10);
            elseif strcmp(obj.Parameter, 'Improved')
                [C, kappa] = obj.PronyCoefficientsImproved(Pressure,MicSpacing,NrModes);
            elseif strcmp(obj.Parameter, 'Improved_2')
                [C, kappa] = obj.PronyCoefficientsImproved_2(Pressure,MicSpacing,NrModes);
            end
            
            error_k = max(abs(1i * MicSpacing * k(1:min([length(k), NrModes, length(kappa)])) - 1i * MicSpacing * kappa(1:min([length(k), NrModes, length(kappa)])))) / max(abs(k));
            error_c = max(abs(A(1:min([length(A), NrModes, length(C)])) - transpose(C(1:min([length(A), NrModes, length(C)]))))) / max(abs(A));
            
            h = 0;
            for ii = 1:min([NrModes, length(kappa)])
                h = h + C(ii) * exp(1i * kappa(ii) * X);
            end
            
            error_h = max(abs(Pressure - h)) / max(abs(Pressure));
            
            disp(error_k)
            disp(error_c)
            disp(error_h)
            
            
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
        function [Amplitude, kappa] = PronyCoefficients(P,dX,NrModes)
            %Input
            %P is the complex microphone pressure which should be a row
            %vectors
            %dX is the microphone seperation
            %NrModes is the amount of modes that have to be decomposed
            M = length(P);
            if NrModes>M/2
                error('The number of modes can not exceed the half of the measurement points')
            end
            %The hankel matrix is constructed, in order to solve for the
            %exponential functions
            H = toeplitz(P(NrModes:(M-1)),fliplr(P(1:(NrModes))));
            h_M = P(NrModes+1:end).';
            %Obtaining the values of the complex exponential functions
            Lambda = roots([1; -pinv(H)*h_M]);
%             Solution = pinv(transpose(H) * H) * transpose(H) * h_M;
%             Lambda = roots([1; -Solution]);
            
            %Construct the vandermonde matrix to obtain the coefficients
            V = Lambda.^(0:M-1).';
            %Obtain the coefficients for each of the exponential
            C = V\P.';
            Amplitude = C;
            %The wavenumbers associated with the coefficients
            kappa = -1i*log(Lambda)/dX;
            
            [C,I] = sort(abs(C),'descend');
            Amplitude=Amplitude(I);
            kappa=kappa(I);
        end
        
        function [Amplitude, kappa] = PronyCoefficientsImproved(P,dX,NrModes)
            %Input
            %P is the complex microphone pressure which should be a row
            %vectors
            %dX is the microphone seperation
            %NrModes is the amount of modes that have to be decomposed
            M = length(P);
            if NrModes>M/2
                error('The number of modes can not exceed the half of the measurement points')
            end
            
%             A = conj(toeplitz(P(NrModes:(M-1)),fliplr(P(1:(NrModes)))));
%             h = conj(P(NrModes+1:end).');
            
            A = toeplitz(P(NrModes:(M-1)),fliplr(P(1:(NrModes))));
            h = P(NrModes+1:end).';
            
            [U, Delta, V] = svd(A);
            
            SingVal = zeros(1, min(size(Delta)));
            
            for ii = 1:min(size(Delta))
                SingVal(ii) = Delta(ii,ii);
            end
            
            not_squared_sum = sum(SingVal);
            
            alpha = 0;
            ii_alpha = 1;
            
            while alpha < 0.999
                alpha = alpha + SingVal(ii_alpha)/not_squared_sum;
                ii_alpha = ii_alpha + 1;
            end
            
            for ii = ii_alpha:length(SingVal)
                Delta(ii, ii) = 0;
            end
            
            A = U*Delta*V.';
           
            Lambda = roots([1; -pinv(A)*h]);
            
            %Construct the vandermonde matrix to obtain the coefficients
            V = Lambda.^(0:M-1).';
            %Obtain the coefficients for each of the exponential
            C = V\P.';
            Amplitude = C;
            %The wavenumbers associated with the coefficients
            kappa = -1i*log(Lambda)/dX;
            
            [C,I] = sort(abs(C),'descend');
            Amplitude=Amplitude(I);
            kappa=kappa(I);
        end
        
        function [Amplitude, kappa] = PronyCoefficientsImproved_2(P,dX,NrModes)
            %Input
            %P is the complex microphone pressure which should be a row
            %vectors
            %dX is the microphone seperation
            %NrModes is the amount of modes that have to be decomposed
            M = length(P);
            if NrModes>M/2
                error('The number of modes can not exceed the half of the measurement points')
            end
            
            A = conj(toeplitz(P(NrModes:(M-1)),fliplr(P(1:(NrModes)))));
            h = conj(P(NrModes+1:end).');
            
%             A = toeplitz(P(NrModes:(M-1)),fliplr(P(1:(NrModes))));
%             h = P(NrModes+1:end).';
            
            [U, Delta, V] = svd(A);
            
%             SingVal = zeros(min(size(Delta)));
            for ii = 1:min(size(Delta))
                SingVal(ii) = Delta(ii,ii);
            end
            
            squared_sum = sum(SingVal.^2);
            
            alpha = 0;
            ii_alpha = 1;
            
            while alpha < 0.8
                alpha = alpha + (SingVal(ii_alpha)^2)/squared_sum;
                ii_alpha = ii_alpha + 1;
            end
            
            for ii = ii_alpha:length(SingVal)
                Delta(ii, ii) = 0;
            end
            
            A_1 = conj(transpose(A)) * A;
            A_2 = A * conj(transpose(A));
            
            [U, EigenValues_u] = eig(A_1);
            [V, EigenValues_v] = eig(A_2);
            
            B = zeros(length(A), 1);
            for ii = 1:min([length(SingVal) M size(U,1) size(V,1)])
                B = B - (1/SingVal(ii)) * (conj(transpose(U(:,ii))) * h) * V(:,ii);
            end
            
            Lambda = roots([1; -B]);
            
            %Construct the vandermonde matrix to obtain the coefficients
            W = Lambda.^(0:M-1).';
            %Obtain the coefficients for each of the exponential
            C = W\P.';
            Amplitude = C;
            %The wavenumbers associated with the coefficients
            kappa = -1i*log(Lambda)/dX;
            
            [C,I] = sort(abs(C),'descend');
            Amplitude=Amplitude(I);
            kappa=kappa(I);
        end
        
        function [Amplitudes, kappa] = MatrixPencil(P, dX, L)
            % L = number of modes that are calculated
            % M = number of modes we look at
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif L < 3
                error('The number of modes must be superior or equal to 3')
            end
            
%             H_s0 = PronyMethod.HankelProny(2*N-L, L, 0, P);
%             H_s1 = PronyMethod.HankelProny(2*N-L, L, 1, P);
            
            H = PronyMethod.HankelProny(2*N-L, L+1, 0, P);
            
            [Q, R, MPi] = qr(H);
            
            M = 1;
            epsilon = 1e-10;
            Rabs = abs(R);
            [size_1, size_2] = size(Rabs);
            for ii = 1:min(size_1, size_2)-1
                if Rabs(M+1, M+1) > epsilon * Rabs(1, 1)
                    if M < L
                        M = M+1;
                    end
                end
            end
            
            S = R * transpose(MPi);
            
%             S_0 = PronyMethod.Matrix_S(S, L, 0);
%             S_1 = PronyMethod.Matrix_S(S, L, 1);
            
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
            
            V = PronyMethod.VandermondeProny(2*N, Z);
            
            C = V \ transpose(P);
            
            Amplitudes = C;
            [C,I] = sort(abs(C),'ascend');
            Amplitudes = Amplitudes(I);
            
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
            kappa = kappa(I);
            
        end
        
        function [Amplitudes, kappa] = MyLittleProny(P, dX, L, epsilon)
            % L = number of modes that are calculated
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            end
            
            H = PronyMethod.HankelProny(2*N-L, L, 0, P);
            
            h = zeros(2*N-L, 1);
            for ii = 1:length(h)
                h(ii) = -P(ii+L);
            end
            
            q = H\h;
            
            Z_0 = roots(q);
            
            V_0 = PronyMethod.VandermondeProny(2*N, Z_0);
            
            C_0 = V_0\transpose(P);
            
            index = [];
            for ii = 1:length(C_0)
                if C_0(ii) > epsilon
                    index = [index, ii];
                end
            end
            
            Z_1 = Z_0(index);
            
            V_1 = PronyMethod.VandermondeProny(2*N, Z_1);
            
            C_1 = V_1\transpose(P);
            
            f_i = log(Z_1);
            
            Amplitudes = C_1;
            [C_1,I] = sort(abs(C_1),'ascend');
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
            
            H = PronyMethod.HankelProny(2*N-L, L+1, 0, P);
            
            [U, D, W] = svd(H);
            
            W = ctranspose(W);
            
            sigma = diag(D);
            a = [];
            
            for ii = 2:length(sigma)
                if sigma(ii) < epsilon * sigma(1)
                    a = [a, ii];
                end
            end
            
            if isempty(a) == 1
                M = length(sigma);
            else
                M = a(1);
            end 
            
            W_rect = W(1:M, 1:L+1);
            
            W_0 = W_rect(1:M, 1:L);
            W_1 = W_rect(1:M, 2:L+1);
            
            F = pinv(transpose(W_0)) * transpose(W_1);
            
            Z = eig(F);
            
            V = PronyMethod.VandermondeProny(2*N, Z);
            
            C = V\transpose(P);
            
            f_i = log(Z);
            
            Amplitudes = C;
            [C,I] = sort(abs(C),'descend');
            Amplitudes = Amplitudes(I);
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
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
        
        function H = HankelProny(M, N, s, P)
            H = zeros(M, N);
            for ii = 1:M
                for jj = 1:N
                    H(ii, jj) = P((ii - 1) + (jj - 1) + s + 1);
                end
            end
        end
        
        function V = VandermondeProny(M, Z)
            for ii = 1:M %M
                    V(ii, :) = Z(:) .^ (ii - 1);
            end
        end
        
        function T = Matrix_T(S, M, L, s)
            T = S(1:M, (1+s):(L+s));
        end
        
        function S_1 = Matrix_S(S_0, L, s)
            S_1 = S_0(1:length(S_0(:,1)), (1+s):(L+s));
        end
        
        function noise = Noise(Length, Amplitude)
            noise_real = randn(1, Length);
            noise_imag = randn(1, Length);
            noise_real = noise_real / max(noise_real);
            noise_imag = noise_imag / max(noise_imag);
            noise = noise_real + 1i * noise_imag;
            noise = Amplitude * noise;
        end
        
    end
    
end
