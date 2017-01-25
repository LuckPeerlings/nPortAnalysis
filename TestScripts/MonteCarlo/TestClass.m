classdef TestClass < matlab.mixin.SetGet
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        A
        z
    end
    
    methods
        function obj = Calculate(obj)
             obj.z(1,:) = sum(obj.A.x - 0.*obj.A.y.^2,1);
%             obj.z(1,2) = sum(obj.A.x(:,2)) + 5*sum(obj.A.y(:,2));
%             obj.z(1,3) = sum(obj.A.x(:,3)) + 5*sum(obj.A.y(:,3));
%             obj.z(1,4) = sum(obj.A.x(:,4)) + 5*sum(obj.A.y(:,4));
        end
    end
    
end

