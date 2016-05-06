classdef MultiVariateAnalysis < handle
    properties
        Input %Structure containing the input property of the the ClassHandle. TO IMPLEMENT: Make this a cell array with a list of properties and their respective values. 
        OutputProperties %Cell array containing the output properties
        ClassHandle
        MethodHandles
        Methods
        Output
        CorrelationInfo
        UVInputList
        UVOutputList
        BaseInput
        BaseOutput
    end
  

    methods 
        function obj = CalculateUncertainty(obj)
            obj = GetListInputUncertainVariables(obj);
            obj = GetListOutputUncertainVariables(obj);
            obj = DetermineSensitivity(obj);
            obj = UnCorrVariance(obj);
            obj = CreateCorrelationInfo(obj);
            obj = CorrVariance(obj);
            obj = SetOutput(obj);
        end
        function obj = SetOutput(obj)
            %Function to save the output to the output properties           
            for ii = 1:length(obj.UVOutputList)                
                if isempty(obj.UVOutputList(ii).StructPos)
                    obj.Output.( obj.OutputProperties{obj.UVOutputList(ii).Pos} ) = obj.UVOutputList(ii).UV;
                else
                    if ~isfield(obj.Output,obj.OutputProperties{obj.UVOutputList(ii).Pos})
                    obj.Output.( obj.OutputProperties{obj.UVOutputList(ii).Pos} ) = [];
                    end
                    obj.Output.( obj.OutputProperties{obj.UVOutputList(ii).Pos} ) = setfield( obj.Output.( obj.OutputProperties{obj.UVOutputList(ii).Pos} ),obj.UVOutputList(ii).StructPos{:},obj.UVOutputList(ii).UV);
                end
            end
        end
        
        obj = GetListInputUncertainVariables(obj); %Fixed
        obj = GetListOutputUncertainVariables(obj);
        Output = UpdateSolution(obj,Input); %Fixed
        obj = CreateCorrelationInfo(obj);
        obj = DetermineSensitivity(obj); %Fixed
        obj = UnCorrVariance(obj); %Fixed
        obj = CorrVariance(obj);    

    end
    methods (Static)
        LocationName = FindObjectInStructure(Input,ObjectName,varargin)
        
        Sensitivity = MultiVariateAnalysis.CalcSensitivity(BaseOutput,PerturbedOutput)
    end


end

%TODO LIST
%2015-11-19 %Merge the two functions findObjectinStructure and
%findValueinStructure from GetListInputUncertainVariables and
%GetListOutputUncertaintVariables
%-------------------------------------------------------------------------