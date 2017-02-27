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
        
        function displayScatMatrix(obj,varargin)
            pars = inputParser;   %Create a parser object
            DEFAULT.Strouhal = [];
            
            % The available arguments and their check function
            addParameter(pars,'Strouhal',DEFAULT.Strouhal)
            
            %Parsing of the arguments
            parse(pars,varargin{:});
            
            %Assigning the parsed arguments to their variables
            Strouhal = pars.Results.Strouhal;
            
            if isempty(Strouhal)
                XCoordinate = obj.ClassHandle.FreqVec;
                XCoordinateLabel = 'Frequency';
            else
                XCoordinate = 2*pi*obj.ClassHandle.FreqVec*Strouhal.L/Strouhal.U;
                XCoordinateLabel = 'Strouhal ka/M';
            end
            
            h_fig_abs = figure;
            FigTitle = annotation(h_fig_abs,'textbox',[0.5-0.1/2 0.9 0.1 0.1]);
            set(FigTitle,'String','Absolute Values')
            set(FigTitle,'HorizontalAlignment','Center')
            set(FigTitle,'LineStyle','None')
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
                    plot(XCoordinate, abs( ScatUV.Value) + 2*sqrt(AlignedVar(1,:)),'k');
                    plot(XCoordinate, abs( ScatUV.Value),'-');
                    plot(XCoordinate, abs( ScatUV.Value) - 2*sqrt(AlignedVar(1,:)),'k');
                    xlabel(XCoordinateLabel)
                    ylabel('[-]')
                end
            end
            h_fig_angle = figure;
            FigTitle = annotation(h_fig_angle,'textbox',[0.5-0.1/2 0.9 0.1 0.1]);
            set(FigTitle,'String','Angles')
            set(FigTitle,'HorizontalAlignment','Center')
            set(FigTitle,'LineStyle','None')
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(jj),num2str(ii)]);
                    subplot(obj.ClassHandle.NrPorts,obj.ClassHandle.NrPorts, (ii - 1)*obj.ClassHandle.NrPorts + jj )
                    hold on
                    plot(XCoordinate, unwrap(angle( ScatUV.Value))*180/pi + 2*sqrt(AlignedVar(4,:))./abs( ScatUV.Value)*180/pi,'k-');
                    plot(XCoordinate, unwrap(angle( ScatUV.Value))*180/pi,'-');
                    plot(XCoordinate, unwrap(angle( ScatUV.Value))*180/pi - 2*sqrt(AlignedVar(4,:))./abs( ScatUV.Value)*180/pi,'k-');
                    xlabel(XCoordinateLabel)
                    ylabel('[Deg]')
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
        
        function exportData_AbsAngle(obj,DataFile)
            fileID = fopen(DataFile,'w');
            
            %Write the headers
            fprintf(fileID,'f');
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    StringCoeff = ['S',num2str(ii),num2str(jj)];
                    fprintf(fileID,['\t',StringCoeff,'_abs',...
                        '\t',StringCoeff,'_angle',...
                        '\t',StringCoeff,'u_abs',...
                        '\t',StringCoeff,'u_angle']);
                end
            end
            fprintf(fileID,'\n');
            
            %Calculate the variance matrix rotated in the direction of
            %the mean vector value.
            %Taken from In-Phase/Quadrature Covariance-Matrix
            %Representation of the Uncertainty of Vectors and Complex
            %Numbers, Dylan F. Willians, C.m. Wand and Uwe Arz
            for ii = 1:obj.ClassHandle.NrPorts
                for jj = 1:obj.ClassHandle.NrPorts
                    ScatUV = obj.Output.ScatNPort.(['S',num2str(ii),num2str(jj)]);
                    if isempty(ScatUV.CorrVar)
                        TotVar = squeeze(sum(ScatUV.Var,1));
                    else
                        TotVar = squeeze(sum(ScatUV.Var,1) + sum(ScatUV.CorrVar,1));
                    end
                    for mm = 1:size(TotVar,1)
                        Theta = angle(ScatUV.Value(1,mm));
                        R = [cos(-Theta), -sin(-Theta); sin(-Theta), cos(-Theta)];
                        AlignedUCMatrix = R*reshape(TotVar(mm,:),2,2)*transp(R);
                        AlignedVar(ii,jj,:,mm) = AlignedUCMatrix(:);
                    end
                end
            end
            %Writing the information to the file
            for ff = 1:length(obj.ClassHandle.FreqVec)
                fprintf(fileID,'%e',obj.ClassHandle.FreqVec(ff));
                for ii = 1:obj.ClassHandle.NrPorts
                    for jj = 1:obj.ClassHandle.NrPorts
                        ScatUV = obj.Output.ScatNPort.(['S',num2str(ii),num2str(jj)]);
                        Angle = unwrap(angle( ScatUV.Value))*180/pi;
                        fprintf(fileID,'\t %e \t %e \t %e \t %e',abs(ScatUV.Value(ff)),Angle(ff), 2*sqrt(AlignedVar(ii,jj,1,ff)),2*sqrt(AlignedVar(ii,jj,4,ff))./abs( ScatUV.Value(ff))*180/pi);
                    end
                end
                fprintf(fileID,'\n');
            end
            fclose(fileID);
        end
        
        function VarianceDistribution_AbsAngle(obj,varargin)
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
                    size(RelativeVar)
                    if ~isempty(varargin)
                        %Save the values of the relative variances in a
                        %large matrix. This will be easier to write the
                        %data to a file then appending each column to the
                        %file.
                        RelativeVar_Export(ii,jj,:,:)= RelativeVar;
                    end
                    area(h1(cc),obj.ClassHandle.FreqVec,RelativeVar.')
                    if ii == 1 && jj == 1
                        legend(UniqueGroupNames)
                        title('Relative Contribution of variance to the magnitude uncertainty (Uncorrelated Errors)')
                    end
                    cc = cc +1;
                end
            end
            if ~isempty(varargin)
                SAVEDIR = varargin{1};
                if ~strcmp(SAVEDIR(end),'\') 
                    SAVEDIR=strcat(SAVEDIR,'\');
                end
                for ii = 1:obj.ClassHandle.NrPorts
                    for jj = 1:obj.ClassHandle.NrPorts
                        %Write each variance distribuion of the scattering
                        %coefficients to a file.
                        fileName = [SAVEDIR,'Variance_Abs_S',num2str(ii),num2str(jj),'.dat'];
                        fileID = fopen(fileName,'w');
                        
                        %Writing header
                        fprintf(fileID,'f \t');
                        for kk = 1:length(UniqueGroupNames)
                            fprintf(fileID,UniqueGroupNames{kk});
                            if kk ~= length(UniqueGroupNames)
                                fprintf(fileID,'\t');
                            else                            
                                fprintf(fileID,'\n');
                            end
                        end
                        %Writing columns to the file 
                        for ff = 1:size(RelativeVar_Export,4)
                            fprintf(fileID,'%e \t',obj.ClassHandle.FreqVec(ff));
                            for kk = 1:size(RelativeVar_Export,3)                                
                                fprintf(fileID,'%e',RelativeVar_Export(ii,jj,kk,ff));
                                if kk ~= size(RelativeVar_Export,3)
                                fprintf(fileID,'\t');
                                else
                                fprintf(fileID,'\n');
                                end
                            end
                        end
                        fclose(fileID);
                    end
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
                    if ~isempty(varargin)
                        %Save the values of the relative variances in a
                        %large matrix. This will be easier to write the
                        %data to a file then appending each column to the
                        %file.
                        RelativeVar_Export(ii,jj,:,:)= RelativeVar;
                    end
                    area(h2(cc), obj.ClassHandle.FreqVec,RelativeVar.')
                    if ii == 1 && jj == 1
                        legend(UniqueGroupNames)
                        title('Relative Contribution of variance to the phase uncertainty (Uncorrelated Errors)')
                    end
                    cc = cc + 1;
                end
            end
            if ~isempty(varargin)
                SAVEDIR = varargin{1};
                if ~strcmp(SAVEDIR(end),'\') 
                    SAVEDIR=strcat(SAVEDIR,'\');
                end
                for ii = 1:obj.ClassHandle.NrPorts
                    for jj = 1:obj.ClassHandle.NrPorts
                        %Write each variance distribuion of the scattering
                        %coefficients to a file.
                        fileName = [SAVEDIR,'Variance_Angle_S',num2str(ii),num2str(jj),'.dat'];
                        fileID = fopen(fileName,'w');
                        
                        %Writing header
                        fprintf(fileID,'f \t');
                        for kk = 1:length(UniqueGroupNames)
                            fprintf(fileID,UniqueGroupNames{kk});
                            if kk ~= length(UniqueGroupNames)
                                fprintf(fileID,'\t');
                            else                            
                                fprintf(fileID,'\n');
                            end
                        end
                        %Writing columns to the file 
                        for ff = 1:size(RelativeVar_Export,4)
                            fprintf(fileID,'%e \t',obj.ClassHandle.FreqVec(ff));
                            for kk = 1:size(RelativeVar_Export,3)                                
                                fprintf(fileID,'%e',RelativeVar_Export(ii,jj,kk,ff));
                                if kk ~= size(RelativeVar_Export,3)
                                fprintf(fileID,'\t');
                                else
                                fprintf(fileID,'\n');
                                end
                            end
                        end
                        fclose(fileID);
                    end
                end
            end
        end
    end
    
end

