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
            for ii = 1:obj.NrPorts
                NrMeas(ii) = length( fields( obj.Input.(['Port',num2str(ii)]) ) );
                if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                    NrMeas(ii) = NrMeas(ii)-1;
                end
            end
            obj.NrMeas = NrMeas;
            
            for ii = 1:obj.NrPorts
                for jj = 1:obj.NrMeas   
                    if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
                       obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = MergeStructs( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]), obj.Input.(['Port',num2str(ii)]).('Constant') );
                    end
                    InputDecomp = obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]);
                    InputDecomp.f = obj.FreqVec;
                    P = NPortAnalysis.WaveDecomposition( InputDecomp);
                    
                    H_R(ii,jj,:) = P.Min(:,:); 
                    H_L(ii,jj,:) = P.Plus(:,:);  
                    
                end
            end
            assignin('base','H_R',H_R)
            assignin('base','H_L',H_L)
            % Save each field of the scattering matrix in the ScatNPort.
            
            for ii = 1:length(obj.FreqVec)
                S(:,:,ii) =  H_L(:,:,ii)\H_R(:,:,ii);
            end
            for ii = 1:size(S,1)
                for jj = 1:size(S,2)
                    obj.ScatNPort.(['S',num2str(jj),num2str(ii)]) = reshape(S(ii,jj,:),1,[]);
                end
            end
            
%             
%             if ~obj.InputChecked
%                 error('The input has not been checked')
%             end
            %Calculate the ingoing and outgoing waves for each port and
            %measurement
            
%             for ii = 1:obj.NrPorts
%                 for jj = 1:obj.NrMeas   
%                     if isfield(obj.Input.(['Port',num2str(ii)]),'Constant')
%                        obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]) = MergeStructs( obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]), obj.Input.(['Port',num2str(ii)]).('Constant') );
%                     end
%                     InputDecomp = obj.Input.(['Port',num2str(ii)]).(['Meas',num2str(jj)]);
%                     InputDecomp.f = obj.FreqVec;
%                     P = NPortAnalysis.WaveDecomposition( InputDecomp);
%                     assignin('base',['P_',num2str(jj)],P)
%                     
% %                     for ff = 1:length(obj.FreqVec)
% %                         if ii == 1
% %                             H_R(:,jj,ff) = P.Min(:,ff); 
% %                             H_L(:,jj,ff) = P.Plus(:,ff);  
% %                         else                            
% %                             H_R(:,jj,ff) = [P.Min(:,ff) H_L(ii-1,jj,ff)]; 
% %                             H_L(:,jj,ff) = [P.Plus(:,ff) H_L(ii-1,jj,ff)];  
% %                         end
% %                     end
%                     %Real Version!!!
%                     H_R(ii,jj,:) = P.Min(1,:); 
%                     H_L(ii,jj,:) = P.Plus(1,:);          
%                 end
%             end

%             %Calculate the scattering matrix as function of frequency   
%             
% %             MaxMode = size(P.Min,1);
% %             if obj.NrMeas == 1
% %                 S = zeros(MaxMode,1,length(obj.FreqVec));
% %                 for ii = 1:length(obj.FreqVec)
% %                     H_R_Stripped = H_R(:,:,ii);
% %                     H_R_Stripped (~all(H_R_Stripped,2), : ) = [];
% %                     H_L_Stripped = H_L(:,:,ii);
% %                     H_L_Stripped (~all(H_L_Stripped,2), : ) = [];
% %                     try
% %                     S_single = H_R_Stripped./H_L_Stripped;
% %                     catch
% %                         S_single = 0;
% %                     end
% %                     S(1:size(S_single,1),1,ii) = S_single;
% %                 end
% %             else
% %                 length(obj.FreqVec)
% %                 S = zeros(MaxMode,MaxMode,length(obj.FreqVec));
% %                 for ii = 1:length(obj.FreqVec)
% %                     H_R_Stripped = H_R(:,:,ii);
% %                     H_R_Stripped (~all(H_R_Stripped,2), : ) = [];
% %                     H_L_Stripped = H_L(:,:,ii);
% %                     H_L_Stripped (~all(H_L_Stripped,2), : ) = [];
% %                     
% %                     S_single = (H_R_Stripped'*H_R_Stripped) \ (H_R_Stripped'*H_L_Stripped)
% %                     S(1:size(S_single,1),1:size(S_single,2),ii) = S_single;
% %                     if obj.NrMeas < size(P.Min,MaxMode)
% %                         warning('The number of independent measurement is not enough to have a full scattering matrix')
% %                         break
% %                     end
% %                 end
% %             end
%             
%             
%             %

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
                    plot(obj.FreqVec, angle( obj.ScatNPort.(['S',num2str(jj),num2str(ii)]) )*180/pi)
                end
            end
        end
            
            
    end
    methods (Static)
        WaveDecomposition_chk(varargin)
        P = WaveDecomposition(varargin)
        WaveNumber_chk(varargin)
        k = WaveNumber(varargin)
        ModalMatrix = RectangularModalMatrix(GasProp,WaveNumberProp,Index)
        ModalMatrix = CircularModalMatrix(GasProp,WaveNumberProp,Index)
        [Gamma_beatty] = getAlphaFromBeatty1950_Rectangular(ny,nz,sh,He,aspectRatio,Pr,gamma)
        [Gamma_beatty] = getAlphaFromBeatty1950_Circular(m,n,sh,He,Pr,gamma,RootsBesselFunction)

    end
   
end
