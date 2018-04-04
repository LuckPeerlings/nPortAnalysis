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
        noise
    end
    
    methods
        %Constructor
        function obj = PronyMethod(Frequency,P,dX,NrModes,Parameter,noise)
            obj.Frequency = Frequency; 
            obj.P = P;
            obj.dX = dX;
            obj.NrModes = NrModes;
            obj.Parameter = Parameter;
            obj.noise = noise;
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
            NrModes = 50;
            MicSpacing = 0.01;
            X = 0:MicSpacing:1-MicSpacing;
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            
            %The wavenumbers are be sorted using the size of the absolute
            %value.
            k = [5+5i,4+4i,2+1i,1]*10; %,2+1i,1
            A = [15+9i,10-8i,9+6i,5-5i]; %,2+2i,1-1i
            
            noise = obj.noise;
            
            %Determine pressure with the above exponentials
            
            Pressure = A(1)*exp(1i*k(1)*X) + A(2)*exp(1i*k(2)*X) + A(3)*exp(1i*k(3)*X) + A(4)*exp(1i*k(4)*X) + (1 + 1i) * noise; %(20 + 20i) * rand(1, length(X)) ;
            if strcmp(obj.Parameter, 'Normal')
%                 [C, kappa] = obj.PronyCoefficients(Pressure,MicSpacing,NrModes);
                [C, kappa] = obj.MatrixPencil(Pressure, MicSpacing, NrModes, NrModes);
            elseif strcmp(obj.Parameter, 'Improved')
                [C, kappa] = obj.PronyCoefficientsImproved(Pressure,MicSpacing,NrModes);
            elseif strcmp(obj.Parameter, 'Improved_2')
                [C, kappa] = obj.PronyCoefficientsImproved_2(Pressure,MicSpacing,NrModes);
            end
            
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
        
        function [Amplitudes, kappa] = MatrixPencil(P, dX, L, M)
            % L = number of modes that are calculated
            % M = number of modes we look at
            
            N = length(P)/2;
            if L > N
                error('The number of modes can not exceed the half of the measurement points')
            elseif L < 3
                error('The number of modes must be superior or equal to 3')
            end
            
            H_s0 = PronyMethod.HankelProny(2*N-L, L, 0, P);
            H_s1 = PronyMethod.HankelProny(2*N-L, L, 1, P);
            
            Z = eig(H_s1, H_s0);
%             [HH_s1, HH_s0, Q, R] = qz(H_s1, H_s0);
%             [AAS,BBS,QS,RS] = ordqz(HH_s1, HH_s0, Q, R,'udo');
%             Z = ordeig(AAS, BBS);
            
            if length(Z) < M
                disp('There are less eigenvalues than modes')
            elseif length(Z) > M
                disp('There are more eigenvalues than modes')
            end
            
%             [ZZ, I] = sort(abs(Z), 'descend');
            f_i = log(Z);
            
%             Z(1:4) = [5+5i,4+4i,2+1i,1]*10;
%             Z = Z(I);
            
            V = PronyMethod.VandermondeProny(2*N, Z);
%             V = Z.^(0:M-1); %(2*N-1));
%             V = transpose(V);
            
            C = V \ transpose(P);
%             C = pinv(V) * transpose(P(1:M+1));
            Amplitudes = C;
            [C,I] = sort(abs(C),'descend');
            Amplitudes = Amplitudes(I);
            
            disp(V * C - transpose(P));
            
            for ii = 1:length(f_i)
                kappa(ii) = -1i * f_i(ii) / dX;
            end
            
%             kappa = kappa(I);
            
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
    end
    
end
