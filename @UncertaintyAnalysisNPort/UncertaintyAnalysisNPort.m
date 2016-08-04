classdef UncertaintyAnalysisNPort < MultiVariateAnalysis
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    
    
    
    properties
        
    end
    
    methods
        function obj = UncertaintyAnalysisNPort()
            obj.OutputProperties = {'ScatNPort'};
            obj.ClassHandle = NPortAnalysis;
            obj.MethodHandles = {'calculateScatteringMatrix'};
        end
        
        function obj = PlotScatMatrix(obj)
            
        end
        
        function displayScatMatrix(obj)
            figure;
            title('Absolute values')
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    if isempty(ScatUV.CorrVar)
                        TotVar = squeeze(sum(ScatUV.Var,1));
                    else
                        TotVar = squeeze(sum(ScatUV.Var,1) + sum(ScatUV.CorrVar,1));
                    end
                    %Calculate the variance matrix rotated in the direction of
                    %the mean vector value.
                    %Taken from In-Phase/Quadrature Covariance-Matrix
                    %Representation of the Uncertainty of Vectors and Complex
                    %Numbers, Dylan F. Willians, C.m. Wand and Uwe Arz
                    for mm = 1:size(TotVar,1)
                        Theta = angle(ScatUV.Value(1,mm));
                        R = [cos(-Theta), -sin(-Theta); sin(-Theta), cos(-Theta)];
                        AlignedUCMatrix = R*reshape(TotVar(mm,:),2,2)*transp(R);
                        AlignedVar(:,mm) = AlignedUCMatrix(:);
                    end
                    subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj );
                    hold on
                    plot(obj.ClassHandle.FreqVec, abs( ScatUV.Value) + 2*sqrt(AlignedVar(1,:)),'k');
                    plot(obj.ClassHandle.FreqVec, abs( ScatUV.Value),'-');
                    plot(obj.ClassHandle.FreqVec, abs( ScatUV.Value) - 2*sqrt(AlignedVar(1,:)),'k');
                    assignin('base','AlignedVar',AlignedVar)
                end
            end
            figure;
            title('Phase Values')
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj )
                    hold on
                    plot(obj.ClassHandle.FreqVec, unwrap(angle( ScatUV.Value))*180/pi + 2*sqrt(AlignedVar(4,:))./abs( ScatUV.Value)*180/pi,'k-');
                    plot(obj.ClassHandle.FreqVec, unwrap(angle( ScatUV.Value))*180/pi,'-');
                    plot(obj.ClassHandle.FreqVec, unwrap(angle( ScatUV.Value))*180/pi - 2*sqrt(AlignedVar(4,:))./abs( ScatUV.Value)*180/pi,'k-');                   
                    
                end
            end
        end
        function VarianceDistribution_RealImag(obj)
            figure
            cc = 1;
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    h1(cc) = subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj );
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    assignin('base','ScatUV',ScatUV)
                    %Get the unique groups
                    UniqueGroupNames = unique([ScatUV.Group{:}]);
                    SummedVar = zeros(length(UniqueGroupNames),size(ScatUV.Var,2),4);
                    for ll = 1:size(ScatUV.Var,1);
                        %Find the index of the position in the unique names
                        I = find(strcmp(UniqueGroupNames,ScatUV.Group{ll}));
                        SummedVar(I,:,:) = SummedVar(I,:,:) + ScatUV.Var(ll,:,:);
                    end
                    RelativeVar = bsxfun(@rdivide,SummedVar(:,:,1),sum(SummedVar(:,:,1),1));
                    area(h1(cc),obj.ClassHandle.FreqVec,RelativeVar.')
                    if ii == 1 && jj == 1
                        legend(UniqueGroupNames)
                        title('Relative Contribution of variance to the magnitude uncertainty (Uncorrelated Errors)')
                    end
                    cc = cc +1;
                end
            end
            
            
            figure  
            cc = 1;
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    h2(cc) = subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj );
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    %Get the unique groups
                    UniqueGroupNames = unique([ScatUV.Group{:}]);
                    SummedVar = zeros(length(UniqueGroupNames),size(ScatUV.Var,2),4);
                    for ll = 1:size(ScatUV.Var,1);
                        %Find the index of the position in the unique names
                        I = find(strcmp(UniqueGroupNames,ScatUV.Group{ll}));
                        SummedVar(I,:,:) = SummedVar(I,:,:) + ScatUV.Var(ll,:,:);
                    end
                    RelativeVar = bsxfun(@rdivide,SummedVar(:,:,4),sum(SummedVar(:,:,4),1));
                    area(h2(cc), obj.ClassHandle.FreqVec,RelativeVar.')
                    if ii == 1 && jj == 1
                        legend(UniqueGroupNames)
                        title('Relative Contribution of variance to the phase uncertainty (Uncorrelated Errors)')
                    end
                    cc = cc + 1;
                end
            end
        end    
        
        function VarianceDistribution_AbsAngle(obj)
            figure
            cc = 1;
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    h1(cc) = subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj );
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    assignin('base','ScatUV',ScatUV)
                    %Get the unique groups
                    UniqueGroupNames = unique([ScatUV.Group{:}]);
                    SummedVar = zeros(length(UniqueGroupNames),size(ScatUV.Var,2),4);
                    for ll = 1:size(ScatUV.Var,1);
                        %Find the index of the position in the unique names
                        I = find(strcmp(UniqueGroupNames,ScatUV.Group{ll}));
                        SummedVar(I,:,:) = SummedVar(I,:,:) + ScatUV.AlignedVar(ll,:,:);
                    end
                    RelativeVar = bsxfun(@rdivide,SummedVar(:,:,1),sum(SummedVar(:,:,1),1));
                    area(h1(cc),obj.ClassHandle.FreqVec,RelativeVar.')
                    if ii == 1 && jj == 1
                        legend(UniqueGroupNames)
                        title('Relative Contribution of variance to the magnitude uncertainty (Uncorrelated Errors)')
                    end
                    cc = cc +1;
                end
            end
            
            
            figure  
            cc = 1;
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    h2(cc) = subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj );
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    %Get the unique groups
                    UniqueGroupNames = unique([ScatUV.Group{:}]);
                    SummedVar = zeros(length(UniqueGroupNames),size(ScatUV.Var,2),4);
                    for ll = 1:size(ScatUV.Var,1);
                        %Find the index of the position in the unique names
                        I = find(strcmp(UniqueGroupNames,ScatUV.Group{ll}));
                        SummedVar(I,:,:) = SummedVar(I,:,:) + ScatUV.AlignedVar(ll,:,:);
                    end
                    RelativeVar = bsxfun(@rdivide,SummedVar(:,:,4),sum(SummedVar(:,:,4),1));
                    area(h2(cc), obj.ClassHandle.FreqVec,RelativeVar.')
                    if ii == 1 && jj == 1
                        legend(UniqueGroupNames)
                        title('Relative Contribution of variance to the phase uncertainty (Uncorrelated Errors)')
                    end
                    cc = cc + 1;
                end
            end
        end    
    end
    
end

