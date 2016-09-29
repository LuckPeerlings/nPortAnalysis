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
            
            %And check if there is enough independent measurements
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
                    obj.ScatNPort.(['S',num2str(jj),num2str(ii)]) = reshape(S(ii,jj,:),1,[]);
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
                    plot(obj.FreqVec, abs( obj.ScatNPort.(['S',num2str(jj),num2str(ii)]) ))
                end
            end
            figure;            
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrPorts
                    subplot(obj.NrPorts,obj.NrPorts, (ii - 1)*obj.NrPorts + jj )
                    plot(obj.FreqVec, unwrap(angle( obj.ScatNPort.(['S',num2str(jj),num2str(ii)]) ))*180/pi)
                end
            end
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
