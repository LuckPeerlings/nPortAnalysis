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
    end
    
    methods
        %Constructor
        function obj = PronyMethod(Frequency,P,dX,NrModes)
            obj.Frequency = Frequency; 
            obj.P = P;
            obj.dX = dX;
            obj.NrModes = NrModes;
        end
        function obj = DetermineImpedance(obj)
            
        end
        function obj = DetermineWaveNumber(obj)
            %Calculating the wavenumber and the amplitude of the modes for
            %each frequency
            for ff = 1:size(obj.P,2)                
                [obj.Amplitude(:,ff),obj.kappa(:,ff)] = obj.PronyCoefficients(obj.P(:,ff).',obj.dX,obj.NrModes);                
            end
        end
        
        function obj = PlotWaveNumber(obj)
            %Function to plot the obtained real and imaginary part of each
            %of the modes and the relative amplitude of the modes compared to the largest mode present. 
            for nn = 1:obj.NrModes
                subplot(obj.NrModes+1,1,nn)
                plot(obj.Frequency,real(obj.Amplitude(nn,:)),'x');
                hold all
                plot(obj.Frequency,imag(obj.Amplitude(nn,:)),'o');
            end
            legend('Real part of wavenumber','Imag. part of wavenumber')  
            
            AX1 = subplot(obj.NrModes+1,1,obj.NrModes+1);            
            for ff = 1:size(obj.Frequency,2)
                plot(obj.Frequency(ff),abs(obj.Amplitude(:,ff))/norm(obj.Amplitude(:,ff)),'x')
                hold on         
                AX1.ColorOrderIndex = 1;
                ylabel('RelAmplitude')
                xlabel('Frequency')
            end
        end
        function TestClass(obj)
            %Function to determine whether the prony method is implemented
            %correctly
            NrModes = 10;
            MicSpacing = 0.01;
            X = 0:MicSpacing:1-MicSpacing;
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            k = [5,3,2,1]*10;
            A = [5,3,2,1];
            %Determine pressure with the above exponentials
            Pressure = A(1)*exp(1i*k(1)*X) + A(2)*exp(1i*k(2)*X) + A(3)*exp(1i*k(3)*X) + A(4)*exp(1i*k(4)*X);
            [C,kappa] = obj.PronyCoefficients(Pressure,MicSpacing,NrModes);
            %Sort the modes by amplitude
            figure
            subplot(3,1,1)
            plot(X,real(Pressure))
            xlabel('Distance')
            ylabel('Real Pressure')
            subplot(3,1,2)
            plot(real(kappa),'x')
            hold all
            plot(k,'o')
            ylabel('Wavenumber')
            xlabel('Mode nr')
            subplot(3,1,3)
            plot(real(A),'x')
            hold all
            plot(real(C),'o')
            ylabel('Amplitude')
            xlabel('Mode nr')
        end
    end
    
    methods(Static)
        function [Amplitude,kappa] = PronyCoefficients(P,dX,NrModes)
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
    end
    
end
