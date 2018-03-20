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
        
        function obj = PlotWaveNumberComplexDomain(obj,ModeNr)
            nn = ModeNr;
            
            [c,I] = sort(obj.Frequency);
            I
            figure
            AX  = axes();
            hold on;
            PronyMethod.ArrowPlotComplexDomain(AX,obj.kappa(nn,I))
            axis square
            
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
            NrModes = 10;
            MicSpacing = 0.01;
            X = 0:MicSpacing:1-MicSpacing;
            %The wavenumbers have to be large enough that there is a
            %variation accross the microphones
            
            %The wavenumbers are be sorted using the size of the absolute
            %value.
            k = [1i,2,2+1i,1]*10;
            A = [15+9i,3-9i,2+2i,1-1i];
            
            %Determine pressure with the above exponentials
            Pressure = A(1)*exp(1i*k(1)*X) + A(2)*exp(1i*k(2)*X) + A(3)*exp(1i*k(3)*X) + A(4)*exp(1i*k(4)*X);
            [C,kappa] = obj.PronyCoefficients(Pressure,MicSpacing,NrModes);
            
            %Sort the modes by amplitude
            figure
            subplot(3,2,1)
            plot(X,real(Pressure));
            xlabel('Distance')
            ylabel('Real Pressure')

            subplot(3,2,2)
            plot(X,imag(Pressure))
            xlabel('Distance')
            ylabel('Imag Pressure')
            
            subplot(3,2,3)
            plot(real(kappa),'x')
            hold all
            plot(real(k),'o')
            ylabel('Real Wavenumber')
            xlabel('Mode nr')
            legend('Input','Educed')
            
            subplot(3,2,4)
            plot(imag(kappa),'x')
            hold all
            plot(imag(k),'o')
            ylabel('Real Wavenumber')
            xlabel('Mode nr')   
            legend('Input','Educed')
                        
            subplot(3,2,5)
            plot(real(A),'x')
            hold all
            plot(real(C),'o')
            ylabel('Real Amplitude')
            xlabel('Mode nr')
            legend('Input','Educed')
            
            subplot(3,2,6)
            plot(imag(A),'x')
            hold all
            plot(imag(C),'o')
            ylabel('Imag Amplitude')
            xlabel('Mode nr')
            legend('Input','Educed')
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
        
        function [AX] = ArrowPlotComplexDomain(AX,Data)
            
            for ii = 1:length(Data)-1
                if ii == 1
                quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'o','filled','MarkerFaceColor','r','MaxHeadSize',0.5)
                else
                                    quiver(real(Data(ii)),imag(Data(ii)),...
                       real(Data(ii+1)-Data(ii)), ...
                       imag(Data(ii+1)-Data(ii)),0,'ko','filled','MaxHeadSize',0.5)
                end                
            end
            plot(real(Data(ii+1)),imag(Data(ii+1)),'rx')
                       
        end
    end
    
end
