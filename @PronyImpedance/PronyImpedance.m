classdef PronyImpedance < PronyMethod
    %UNTITLED Summary of this class goes here
    properties (Access = public)
        % Settable properties
        FlowVelocity
        Temperature
        Height
        HeightY
    end    
    properties (Access = protected)
        Density   
        SpeedOfSound
        MachNumber
        FreeFieldWaveNumber
        WaveNumberY %Input
        WaveNumberX
        %Input
        WaveNumberZ 
    end    
    properties (GetAccess = public, SetAccess = protected)
        Impedance
    end
    
    methods
        
        function obj = PronyImpedance(  Frequency,P,MicPositions,MicEqPositions,NrModes,Epsilon,...
                                        FlowVelocity, Temperature, Height)
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
            obj.HeightY = Height;            
            obj.FlowVelocity = FlowVelocity;
            obj.Temperature = Temperature;
        end
        function obj = CalculateImpedance(obj)
            obj.EquispaceData();
            obj.DetermineWaveNumber();
            Properties = AirProperties('t',obj.Temperature);
            obj.Density = Properties.Density;
            obj.SpeedOfSound = Properties.SpeedOfSound;
            obj.MachNumber = obj.FlowVelocity/obj.SpeedOfSound;
           
            obj.CalculateWaveNumbers()
            obj.IngerdMeyersBoundaryCondition()
%             obj.PlotImpedance();
        end
        
        
        function obj = CalculateWaveNumbers(obj,n)
            if nargin == 1
                n = 0;
            end
            obj.FreeFieldWaveNumber = 2*pi*obj.Frequency./obj.SpeedOfSound;
            obj.WaveNumberX = 0;
            obj.WaveNumberZ = obj.WaveNumber + n*2*pi/obj.MicSpacing;
            obj.WaveNumberY = sqrt((obj.FreeFieldWaveNumber-obj.MachNumber.*obj.WaveNumberZ).^2-obj.WaveNumberX.^2-obj.WaveNumberZ.^2);                    
        end
        
        function obj = IngerdMeyersBoundaryCondition(obj)
            obj.Impedance =  (obj.WaveNumberY.*tan(obj.WaveNumberY.*obj.HeightY)).^-1 .* 1i./obj.FreeFieldWaveNumber.*(obj.FreeFieldWaveNumber-obj.MachNumber.*obj.WaveNumberZ).^2; 
            
        end        
        
        function obj = PlotImpedance(obj)
            figure
            set(gcf,'Position',[678 732 987 366]);
            AX = subplot(1,2,1);
            hold on;
            PronyMethod.ArrowPlotComplexDomain(AX,obj.WaveNumberZ);
            axis square;
            grid on
            ylabel('Imag. part of wavenumber')  
            xlabel('Real part of wavenumber')
            
            subplot(1,2,2);
            plot(obj.Frequency, real(obj.Impedance));
            grid on
            hold all
            plot(obj.Frequency, imag(obj.Impedance));
            xlabel('Frequency [Hz]')
            ylabel('Relative Impedance [-]')
            legend('Real part','Imaginary part')
        end
    end
    
end

