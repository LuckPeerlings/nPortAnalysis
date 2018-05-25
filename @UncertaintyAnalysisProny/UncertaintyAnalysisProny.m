classdef UncertaintyAnalysisProny < MultiVariateAnalysis
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = CalculateUncertainty(obj)
            obj.GetListInputUncertainVariables;
            obj.GetListOutputUncertainVariables;
            obj.DetermineSensitivity;
            obj.UnCorrVariance;
            obj.CreateCorrelationInfo;
            obj.CorrVariance;
            obj.SetOutput;
        end
        
        function obj = UncertaintyAnalysisProny(obj)        

        end
        
        function obj = PlotWavenumber(obj)
        
        end
        
        
        function obj = PlotWaveNumberComplexDomain(obj,ModeNr)
            if nargin == 1
                nn = 1;
            else
                nn = ModeNr;
            end                    
            [~,I] = sort(obj.ClassHandle.Frequency);
            obj.Output.WaveNumber
            
            figure;
            AX  = axes();
            hold on;
            PronyMethod.ArrowPlotComplexDomain(AX,obj.Output.WaveNumber.Value(I))
            axis square;
            obj.Output.WaveNumber.plotUncertaintyEllipse(AX,2,[0.5,0,0]);
            
            ylabel('Imag. part of wavenumber')  
            xlabel('Real part of wavenumber')
        end
        
        
        
    end
    
end


