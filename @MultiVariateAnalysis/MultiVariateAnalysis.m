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
        ClassHandleProtected
    end

    methods 
        function obj = CalculateUncertainty(obj)
            %obj = CreateInput;
            obj = GetListInputUncertainVariables(obj);
            obj = GetListOutputUncertainVariables(obj);
            obj = DetermineSensitivity(obj);
            obj = UnCorrVariance(obj);
            obj = CreateCorrelationInfo(obj);
            obj = CorrVariance(obj);
            obj = SetOutput(obj);
        end
        
        function obj = CreateInput(obj)
            %Get the class meta data
            mco = metaclass(obj.ClassHandle);
            %Go over the propertylist and find the inputs and add the value and name to the Input property
            %Go over the propertylist and find the outpus and add these to
            %the OutputProperties
            PropertiesList = mco.PropertyList;
            nn = 1;
            mm = 1;
            for ii = 1:length(PropertiesList)
                if strcmp(PropertiesList(ii).SetAccess,'public')
                    obj.Input{nn}{1} = PropertiesList(ii).Name;
                    UV = obj.ClassHandle.(PropertiesList(ii).Name);
                    obj.Input{nn}{2} = UV;
                    nn = nn +1;
                end                
                if strcmp(PropertiesList(ii).GetAccess,'public') && strcmp(PropertiesList(ii).SetAccess,'protected')
                    obj.OutputProperties{mm} = PropertiesList(ii).Name;
                    mm = mm +1;
                end
            end
            obj.ClassHandleProtected = copy(obj.ClassHandle);
        end
        
        function obj = SetOutput(obj)
            %Function to save the output to the output properties           
            for ii = 1:length(obj.UVOutputList)   
                %Calculate the total uncertainty for each of the uncertain
                %variables.
                obj.UVOutputList(ii).UV.calculateTotalUncertainty;
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