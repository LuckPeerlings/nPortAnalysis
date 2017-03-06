classdef NPortAnalysis  < matlab.mixin.SetGet
    properties
        FreqVec
        Input
    end
    
    properties
        NrPorts = 2;
        NrMeas = 2;
        InputChecked = true;
        ScatNPort
        DissipationNPort
    end
    methods
        function obj = checkInput(obj)
            % Function to check if the provided input is correct.
            % First the obtain the number of ports and independent
            % measurements
            NrPorts = length(fields(obj.Input));
            for ii = 1: NrPorts
                NrMeas(ii) = length( fields( obj.Input.(['Port',num2str(ii)]) ) );
                if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                    NrMeas(ii) = NrMeas(ii)-1;
                end
            end
            
            %Check that the number of measurements is equal for all ports
            NrMeas = unique(NrMeas);
            if length(NrMeas) ~= 1
                error('The number of measurements per port are not equal for all the ports')
            end
            
            %And check if there are enough independent measurements
            if  NrPorts > NrMeas
                error('The number of independent measurements is smaller than the number of ports')
            end
            
            %Parse the input through the various check functions
            for ii = 1:length(NrPorts)
                
                %Check if the input structure has an constant field. This
                %indicates that that information is constant for all the
                %measurements.
                for jj = 1:length(NrMeas)
                    %If the constant field is present, save the constant
                    %fields in the input data.
                    if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                        obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = MergeStructs( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]), obj.Input.(['Port',num2str(ii)]).Constant );
                        
                    end
                    GasProperties_chk(  obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).GasProp)
                    obj.WaveNumber_chk( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).WaveNumberProp, ...
                        'GasProp',obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).GasProp,...
                        'f',obj.FreqVec);
                    obj.WaveDecomposition_chk( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) ,...
                        'GasProp',obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).GasProp,...
                        'f',obj.FreqVec)
                end
            end
            %Finally save the number of ports, measurements and set the
            %InputCheckedFlag
            obj.NrPorts  = NrPorts;
            obj.NrMeas = NrMeas;
            obj.InputChecked = true;
        end       
        function obj = EnergyDissipationPerPort(obj,varargin)
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrMeas
                    if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                        InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = MergeStructs( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]), obj.Input.(['Port',num2str(ii)]).('Constant') );
                    else
                        InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]);
                    end
                end
            end
            %Obtain the mean flow velocity for the ports.
            MachNumberPort = zeros(1,obj.NrPorts);
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrMeas
                    Properties = GasProperties(InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).GasProp);
                    %Take the average of the the measured speeds
                    if isfield(InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).WaveNumberProp,'U')
                        MachNumberPort(ii) = MachNumberPort(ii) + InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).WaveNumberProp.U/mean(Properties.SpeedOfSound)/obj.NrMeas;
                    else
                        warning('No velocity data found');
                    end
                    %Add velocity correction if necessary
                    if isfield(InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]),'Corr')
                        if isfield(InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).Corr,'U')
                            MachNumberPort(ii) = MachNumberPort(ii) + InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).Corr.U/mean(Properties.SpeedOfSound)/obj.NrMeas;
                        end
                    end
                end
                Radius(ii) = InputDecomp.(['Port',num2str(ii)]).Meas1.WaveNumberProp.Model.r;
            end
            %Create the matrices to convert the propagating pressure wave to the propagating
            Mach_Plus = zeros(obj.NrPorts);
            Mach_Min = zeros(obj.NrPorts);
            for ii = 1:obj.NrPorts
                Mach_Plus(ii,ii) = (1-MachNumberPort(ii))*sqrt(pi*Radius(ii)^2);
                Mach_Min(ii,ii) =  (1+MachNumberPort(ii))*sqrt(pi*Radius(ii)^2);
            end
            %Reshape the scattering matrix and calculate the scattering matrix
            %w.r.t to power.
            for ff = 1:length(obj.FreqVec)
                for ii = 1:obj.NrPorts
                    for jj = 1:obj.NrPorts
                        S(ii,jj) = obj.ScatNPort.(['S',num2str(ii),num2str(jj)])(ff);
                    end
                end                           
                Sq = Mach_Min*S*inv(Mach_Plus);
                P_dis(ff,:,:) = eye(2)-transp(conj(Sq))*Sq;
            end
                       
            for ii = 1:obj.NrPorts
                obj.DissipationNPort.(['Port',num2str(ii)]) = squeeze(P_dis(:,ii,ii)).';
            end
            
            %Show the power dissipated per port
            if ~isempty(varargin) && strcmp(varargin{1},'plot')
            figure;
            plot(obj.FreqVec,P_dis(:,1,1),'.','MarkerSize',10)
            hold on        
            plot(obj.FreqVec,P_dis(:,2,2),'r.','MarkerSize',10)
            legend('Port1','Port2')
            end
            if ~isempty(varargin) && length(varargin) == 2
                fileName = varargin{2};
                
                %Write each variance distribuion of the scattering
                %coefficients to a file.
                fileID = fopen(fileName,'w');
                
                %Writing header
                fprintf(fileID,'f \t');
                for ii = 1:obj.NrPorts
                    fprintf(fileID,['Port',num2str(ii),'\t']);
                end
                fprintf(fileID,'\n');
                %Writing columns to the file
                for ff = 1:length(obj.FreqVec)
                    fprintf(fileID,'%e \t',obj.FreqVec(ff));
                    for ii = 1:obj.NrPorts
                        fprintf(fileID,'%e',P_dis(ff,ii,ii));
                        if ii ~= obj.NrPorts
                            fprintf(fileID,'\t');
                        else
                            fprintf(fileID,'\n');
                        end
                    end
                end
                fclose(fileID);
            end
        end       
        function obj = calculateDissipation(obj)
            %Obtain the mean value of the flow speed for each port.
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrMeas
                    if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                        InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = MergeStructs( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]), obj.Input.(['Port',num2str(ii)]).('Constant') );
                    else
                        InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]);
                    end
                end
            end
            %Obtain the mean flow velocity for the ports.
            MachNumberPort = zeros(1,obj.NrPorts);
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrMeas
                    Properties = GasProperties(InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).GasProp);
                    %Take the average of the the measured speeds
                    MachNumberPort(ii) = MachNumberPort(ii) + InputDecomp.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).WaveNumberProp.U/mean(Properties.SpeedOfSound)/obj.NrMeas;
                end
                Radius(ii) = InputDecomp.(['Port',num2str(ii)]).Meas1.WaveNumberProp.Model.r;
            end
            
            %Create the matrices to convert the propagating pressure wave to the propagating
            Mach_Plus = zeros(obj.NrPorts);
            Mach_Min = zeros(obj.NrPorts);
            for ii = 1:obj.NrPorts
                Mach_Plus(ii,ii) = (1-MachNumberPort(ii))*sqrt(pi*Radius(ii)^2);
                Mach_Min(ii,ii) = (1+MachNumberPort(ii))*sqrt(pi*Radius(ii)^2);
            end
            
            %Reshape the scattering matrix and calculate the scattering matrix
            %w.r.t to power.
            for ff = 1:length(obj.FreqVec)
                for ii = 1:obj.NrPorts
                    for jj = 1:obj.NrPorts
                        S(ii,jj) = obj.ScatNPort.(['S',num2str(ii),num2str(jj)])(ff);
                    end
                end
                S_Power = Mach_Min*S*inv(Mach_Plus);
                lambda(:,ff) = eig(S_Power'*S_Power);
                
            end
            figure; plot((1-lambda).')
        end
        function obj = calculateScatteringMatrix(obj)
            obj.NrPorts = length(fields(obj.Input));
            %Check how many measurements there are, if there is a field
            %called constant, do not take this as a measurement
            for ii = 1:obj.NrPorts
                NrMeas(ii) = length( fields( obj.Input.(['Port',num2str(ii)]) ) );
                if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                    NrMeas(ii) = NrMeas(ii)-1;
                end
            end
            assignin('base','Input',obj.Input)
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrMeas
                    if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                        InputDecomp = MergeStructs( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]), obj.Input.(['Port',num2str(ii)]).('Constant') );
                    else
                        InputDecomp = obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]);
                    end
                    InputDecomp.f = obj.FreqVec;
                    [P,Correction] = NPortAnalysis.WaveDecomposition( InputDecomp);
                    
                    
                    if ~isempty(Correction)
                        %If the wave decomposition has an optimization
                        %procedure, save the correction to the input of the
                        %object.
                        obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]).Corr = Correction;
                    end
                    H_R(ii,jj,:) = P.Plus(:,:);
                    H_L(ii,jj,:) = P.Min(:,:);
                    
                end
            end
            % Save each field of the scattering matrix in the ScatNPort.
            
            for ii = 1:length(obj.FreqVec)
                S(:,:,ii) =  H_L(:,:,ii)/H_R(:,:,ii);
                
                [msgstr, msgid] = lastwarn;
                if strcmp(msgid,'MATLAB:illConditionedMatrix')
                    fprintf('Ill conditioned matrix at frequency %f \n',obj.FreqVec(ii))
                    lastwarn('')
                end
            end
            for ii = 1:size(S,1)
                for jj = 1:size(S,2)
                    obj.ScatNPort.(['S',num2str(ii),num2str(jj)]) = reshape(S(ii,jj,:),1,[]);
                end
            end
        end       
        function displayScatMatrix(obj)
            if ~obj.InputChecked
                error('The input has not been checked')
            end
            
            figure;
            
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrPorts
                    subplot(obj.NrPorts,obj.NrPorts, (ii - 1)*obj.NrPorts + jj )
                    plot(obj.FreqVec, abs( obj.ScatNPort.(['S',num2str(ii),num2str(jj)]) ))
                end
            end
            figure;
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrPorts
                    subplot(obj.NrPorts,obj.NrPorts, (ii - 1)*obj.NrPorts + jj )
                    plot(obj.FreqVec, unwrap(angle( obj.ScatNPort.(['S',num2str(ii),num2str(jj)]) ))*180/pi)
                end
            end
        end      
        function displayTransferMatrix(obj)
            %Calculate the scattering transfer matrix
            for ff = 1:length(obj.FreqVec)
                Tm(1,1) = -1/obj.ScatNPort.S12(ff)*(obj.ScatNPort.S22(ff)*obj.ScatNPort.S11(ff)-obj.ScatNPort.S12(ff)*obj.ScatNPort.S21(ff));
                Tm(1,2) =  1/obj.ScatNPort.S12(ff)*obj.ScatNPort.S22(ff);
                Tm(2,1) = -1/obj.ScatNPort.S12(ff)*obj.ScatNPort.S11(ff);
                Tm(2,2) =  1/obj.ScatNPort.S12(ff);

                ConvMat_r = [1 1; 1 -1];
                T(:,:,ff) = ConvMat_r*Tm*inv(ConvMat_r);
                assignin('base','T',T)
            end
             
%             figure; plot(obj.FreqVec, real(Z_1)); hold on; plot(obj.FreqVec, real(Z_2),'-.')
%             figure; plot(obj.FreqVec, imag(Z_1)); hold on; plot(obj.FreqVec, imag(Z_2),'-.')
        end
        function exportData_AbsAngle(obj,DataFile)
            fileID = fopen(DataFile,'w');
            
            %Write the headers
            fprintf(fileID,'f');
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrPorts
                    StringCoeff = ['S',num2str(ii),num2str(jj)];
                    fprintf(fileID,['\t',StringCoeff,'_abs',...
                        '\t',StringCoeff,'_angle']);
                end
            end
            fprintf(fileID,'\n');
            
            
            %Writing the information to the file
            for ff = 1:length(obj.FreqVec)
                fprintf(fileID,'%e',obj.FreqVec(ff));
                for ii = 1:obj.NrPorts
                    for jj = 1:obj.NrPorts
                        ScatElement = obj.ScatNPort.(['S',num2str(ii),num2str(jj)]);
                        Angle = unwrap(angle(ScatElement))*180/pi;
                        fprintf(fileID,'\t %e \t %e \t %e \t %e',abs(ScatElement(ff)),Angle(ff));
                    end
                end
                fprintf(fileID,'\n');
            end
            fclose(fileID);
        end        
    end
    methods (Static)
        WaveDecomposition_chk(varargin)
        [P,Correction] = WaveDecomposition(varargin)
        WaveNumber_chk(varargin)
        k = WaveNumber(varargin)
        ModalMatrix = RectangularModalMatrix(GasProp,WaveNumberProp,Index)
        ModalMatrix = CircularModalMatrix(GasProp,WaveNumberProp,Index)
        [Gamma_beatty] = getAlphaFromBeatty1950_Rectangular(ny,nz,sh,He,aspectRatio,Pr,gamma)
        [Gamma_beatty] = getAlphaFromBeatty1950_Circular(m,n,sh,He,Pr,gamma,RootsBesselFunction)
        
    end
    
end
