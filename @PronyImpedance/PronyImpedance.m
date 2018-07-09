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
        WaveNumberY 
        WaveNumberX
        WaveNumberZ 
    end
    properties (GetAccess = public, SetAccess = protected)
        Impedance
    end
    properties
        UC_Attributes_PronyImpedance = {'FlowVelocity','Input';}
    end
            
    methods
        
        function obj = PronyImpedance(  Frequency,P,MicPositions,MicEqPositions,NrModes,Epsilon,...
                                        FlowVelocity, Temperature, Height)
            obj.Frequency = Frequency; 
            obj.P = P;
            obj.MicPositions = MicPositions;
            obj.MicEqPositions = MicEqPositions;
            obj.NrModes = NrModes;
            obj.Epsilon = Epsilon; 
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
%           obj.PlotImpedance();
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
            set(gcf,'Position',[100 100 987 366]);
            AX = subplot(1,2,1);
            hold on;
            PronyMethod.PlotComplexDomain(AX,obj.WaveNumberZ);
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
        
        
        function obj = ExportToExcel(obj,FileName,Comments)
            %Save Information  on the measurement 
            A{1,1} = 'Comments: ';
            for ii = 1:length(Comments)
                A{ii,2} = Comments{ii};
            end
            
            B = {'','';...
                 'Measurement Conditions','';...                 
                 'Mean Flow Velocity [m/s] ', abs(obj.FlowVelocity);...
                 'Speed of Sound [m/s]', obj.SpeedOfSound;...
                 'Density [kg/m^3]', obj.Density;...
                 '','';...
                 'Algorithm Details','';...
                 'Method ', obj.Method;...
                 'Epsilon', obj.Epsilon;...
                 'Boundary Condition', 'Ingard-Myers';...
                 '','';...
                 };
                      
            for ii = 1:size(B,1)
                    A{ii+length(Comments),1} = B{ii,1};
                    A{ii+length(Comments),2} = B{ii,2};
                
            end
            for ii = 1:size(A,1)
                for jj = 3:5 
                    A{ii,jj} = ''; 
                end
            end
            A(end + 1,:) = {'Frequency [Hz]', 'Real Part Wavenumber [m^-1]', 'Imaginary Part Wavenumber [m^-1]', 'Real Part Normalized Impedance [-]', 'Imaginary Part Normalized Impedance [-]'};
            N = size(A,1);
            for ii = 1:length(obj.Frequency)
                A{ii+N,1} = obj.Frequency(ii);
                A{ii+N,2} = real(obj.WaveNumber(ii));
                A{ii+N,3} = imag(obj.WaveNumber(ii));
                A{ii+N,4} = real(obj.Impedance(ii));
                A{ii+N,5} = imag(obj.Impedance(ii));
            end
            xlswrite(FileName,A)
        end
        
        function obj = ExportToPGF(obj,FileName,Comments)
            Header = {'Frequency', 'RealWavenumber', 'ImaginaryWavenumber', 'RealImpedance', 'ImaginaryImpedance'};
            Data =  [obj.Frequency; real(obj.WaveNumber); imag(obj.WaveNumber); real(obj.Impedance); imag(obj.Impedance)];
            WriteToTextFile_PGF(FileName, Header, Data.', Comments)
        end
    end
    
end

