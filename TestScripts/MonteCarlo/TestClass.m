classdef TestClass < matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        z
    end
    
    methods
        function obj = Calculate(obj)
             obj.z = (obj.A-1-1i);
        end
    end
    
end

