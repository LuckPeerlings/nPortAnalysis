close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Script to use the VXI system to measure and calculate the microphone
%   calibration data. The scrip uses the me4x driver.
%   
%   Luck Peerlings 2014-07-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CHANNEL SETUP
% The microphones should be attached to the first channels and the
% reference to the next channel. So for 2 microphones, chan1 mic1, chan2
% mic2, chan3 ref
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Welcome to the ultimate Calibrator™ 3000!')

[DATA_FILE,DATA_DIR] = uiputfile('.mat','Give the file to save the calibration data');
Data_File = [DATA_DIR,DATA_FILE];


fprintf('To wich channels are the microphones attached? \n')
fprintf('When all channels have been given, do not enter a number\n')

ii = 1;
while  ii<40
    fprintf('Microphone %i \t \t \t',ii)
    ChanInput =                input('Chan: ?                       ');
    if isempty(ChanInput)
        break;
    end
    INPUT_CHAN(ii) = ChanInput;
    ii =  ii+1;
end

REF_SOURCE_CHAN=       input('To which channel is the reference source attached?           '); 
REFERENCE_MIC =        input('To which channel is the reference microphone attached?       ');
SOURCE_CHAN =          input('What is the source channel?                                  ');
ICP =                  input('Are they ICP Microphones [Y/N]?                              ','s');
SENSITIVITY_NEXUS =    input('What is the sensitivity of the NEXUS system [mV/Pa]?         ');
CALIBRATION_PRESSURE = input('What is the needed calibration Pressure [Pa]?                ');
AVERAGES =             input('How many averages are needed [-]?                            ');




INPUT_CHAN = [INPUT_CHAN REF_SOURCE_CHAN];
disp('Which Microphone (ID) is attached to channel:')
for ii = 1:length(INPUT_CHAN)-1
    INPUT_NAMES{ii} = input(['#',num2str(INPUT_CHAN(ii)),' :?    '],'s');
end
INPUT_NAMES{ii+1} = 'REF';

%Start and stop frequency for the calibration data and the number of points
%within this interval
STARTFREQ = 50;
STOPFREQ = 4500;
NR_FREQ = 449;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Starting the measurement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% 1.1 Opening the device
h = actxserver('me4x.MeasEngine');
invoke(h,'Open');
set(h,'MeasType','Sine');


%%%%% 1.2 Initialize inputs    
%Disable all Inputs
set(get(h,'Input'),'Enable',0,0);
% Enable Inputs and set AutoRange.
% InputChannelNumbers
Input_Chan = INPUT_CHAN;
Input_Names = INPUT_NAMES;
for i = 1:length(Input_Chan)
    set(get(h,'Input'),'Enable',Input_Chan(i),1);
    set(get(h,'Input'),'Name',Input_Chan(i),Input_Names{i});
    set(get(h,'Input'),'Range',Input_Chan(i),2);
end
if strcmp(ICP,'Y')
    for i = 1:length(Input_Chan)-1
    set(get(h,'Input'),'Function',Input_Chan(i),'IEPE');
    end
end

%%%%% 1.3 Initialize Outputs
set(get(h,'Source'),'Enable',0,0);
Output_Chan = SOURCE_CHAN;
Output_Input_Chan = REF_SOURCE_CHAN;
%Start amplitude of the calibration.
Output_Amplitude = 0.4; 
Source_Name = {'S1'};
for i = 1:length(Output_Chan)
    set(get(h,'Source'),'Enable',Output_Chan(i),1);
    set(get(h,'Source'),'Name',Output_Chan(i),Source_Name{i});
    set(get(h,'Source'),'Amplitude',Output_Chan(i),Output_Amplitude(i));
end

%%%%% 1.4 Swept sine properties
set(get(h,'Sine'),'StartFreq',STARTFREQ);
set(get(h,'Sine'),'StopFreq',STOPFREQ);
set(get(h,'Sine'),'SweepMode','Linear');
set(get(h,'Sine'),'SweepDirection','Down');
set(get(h,'Sine'),'SweepFreqlines',NR_FREQ); %Max 32767

% Set if the MeasEngine should choose the frequency resolution, dependent
% on the gradient of the response
set(get(h,'Sine'),'AutoResolution',0);

% The integration filter’s effective bandwidth is inversely proportional to
% the integration time (band width = 1 / integration time.)     Increasing
% integrate time effectively narrows the bandwidth at each measurement
% point.  The result is greater harmonic rejection and increased
% signal-to-noise ratios but longer measurement times.

 set(get(h,'Sine'),'IntegrationMode','Cycle');
 set(get(h,'Sine'),'IntegrationTime',50); 
 set(get(h,'Sine'),'SettlingMode','Cycle');
 set(get(h,'Sine'),'SettlingTime',5);

% Properties to use a reference channel.
Ref_Channel = 1;
Ref_Level = 20*log10(SENSITIVITY_NEXUS/1000*CALIBRATION_PRESSURE); %in [dbV] -> 0 = 1V
Ref_MaxLevel = 20*log10(10.0); %in [dbV] -> 0 = 1V
Ref_Tolerance = 20*log10(0.01); %in [dbV] -> 0 = 1V

set(get(h,'Sine'),'SweepDirection','Down');
set(get(h,'Sine'),'AutoLevel',1);
set(get(h,'Sine'),'AutoLevelRefChannel',Ref_Channel);
set(get(h,'Sine'),'AutoLevelRefLevel',Ref_Level);
set(get(h,'Sine'),'AutoLevelMaxLevel',Ref_MaxLevel);
set(get(h,'Sine'),'AutoLevelTolerance',Ref_Tolerance);

%%%%% 2.0 Begin of the measurement loop
NrRepetitions = AVERAGES;

% Create figures before the measurement to increase speed.
cc = hsv(12);
figure
    AX1 = subplot(2,1,1); hold; 
    ylabel(AX1, 'Magnitude',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    xlabel(AX1, 'Frequency',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    AX2 = subplot(2,1,2); hold;
    ylabel(AX2, 'Phase [Rad]',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    xlabel(AX2, 'Frequency',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    set([AX1,AX2], ...
            'Box'           , 'off'     , ...
            'TickDir'       , 'out'     , ...                
            'XScale'        ,'linear'      ,...      
            'Fontsize'      ,16         ,...
            'FontWeight'    ,'Bold',...
            'FontName','Times',...
            'XMinorTick'    ,'on'       ,...
            'YMinorTick'    ,'on'       ,...
            'XColor'      , [0 0 0], ...
            'YColor'      , [0 0 0], ...
            'LineWidth'     ,0.5,...
            'Xgrid'         ,'off',...
            'Ygrid'         ,'off',... 
            'Box'           ,'on'...
        );
    set(gcf, 'color', 'white');

%%%%% 2.1 Measurement
disp('Starting calibration')
for k=1:NrRepetitions
    for i = 1:length(Output_Chan)
        set(get(h,'Source'),'Amplitude',Output_Chan(i),Output_Amplitude(i));
    end 
    invoke(h,'Start');
    disp(['	Repetition ',num2str(k)])
    while 1  
       % Check if new data is available.       
        if get(h,'DataAvailable');
          x = get(h,'xData', 0);
          if length(x) > 0
            clear X Y;  
            for i = 1:length(Input_Chan)  
                X(:,i) = get(h,'xData', 0)';
                Y(:,i) = get(h,'yReal',Input_Chan(i))' + 1i * get(h,'yImag',Input_Chan(i))';
            end        
            if exist('hPlotData1','var')
            delete(hPlotData1)
            delete(hPlotData2)
            end                            
            % Plot the data
            for i = 1:length(Input_Chan)-1  
                frf = Y(:,i)./Y(:,end);    
                hPlotData1(:,i) = plot(AX1,X(:,i), abs(frf), 'color',cc(i,:));
                hPlotData2(:,i) = plot(AX2,X(:,i), angle(frf), 'color',cc(i,:));
            end      
            drawnow;
          end       
          % Notify that I am done with this data.
          invoke(h,'ReadDone');
        end
        if strcmp(get(h,'MeasState'),'ME_READY_FOR_INIT_STATE') && ~get(h,'DataAvailable')
            break
        end
    end
    invoke(h,'Stop');
    
    MeasurementInstance = ['Nr',num2str(k)];
    Measurement.(MeasurementInstance).freq = X;
    Measurement.(MeasurementInstance).Y = Y;
    
end
disp('Calibration has finished')
Measurement.NrMics = length(INPUT_CHAN)-1;
Measurement.Averages = AVERAGES;
Measurement.SensitvityNexus = SENSITIVITY_NEXUS;
Measurement.CalibrationPressure = CALIBRATION_PRESSURE;
Measurement.ReferenceChannel = REFERENCE_MIC;
Measurement.ChannelNames = INPUT_NAMES;
Measurement.Date = datestr(clock);
Measurement.ReferenceMic = INPUT_NAMES{REFERENCE_MIC==INPUT_CHAN};

% Disconnect the VXI hardware.
invoke(h,'Close');

% Cleanup COM objects.
release(h);

% for k=1:NrRepetitions
%     MeasurementInstance = ['Nr',num2str(k)];
%     if sum(Measurement.(MeasurementInstance).freq == 50 ) == 1;
%         warning('Measurement data contains a point at 50Hz, removing that point due to electrical noise')
%         Measurement.(MeasurementInstance).Y(Measurement.(MeasurementInstance).freq == 50,:) = [];
%         Measurement.(MeasurementInstance).f(Measurement.(MeasurementInstance).freq == 50) = [];
%     end
%     if sum(Measurement.(MeasurementInstance).freq == 150 ) == 1;
%         warning('Measurement data contains a point at 150Hz, removing that point due to electrical noise')
%         Measurement.(MeasurementInstance).Y(Measurement.(MeasurementInstance).freq == 150,:) = [];
%         Measurement.(MeasurementInstance).f(Measurement.(MeasurementInstance).freq == 150) = [];
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculating the Calibration data
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
    AX3 = subplot(2,1,1); hold; 
    ylabel(AX3, 'Magnitude',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    xlabel(AX3, 'Frequency',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    AX4 = subplot(2,1,2); hold;
    ylabel(AX4, 'Phase [Rad]',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    xlabel(AX4, 'Frequency',...
           'Fontsize',20, ...
           'Fontweight', 'bold'...
           );
    set([AX3,AX4], ...
            'Box'           , 'off'     , ...
            'TickDir'       , 'out'     , ...                
            'XScale'        ,'linear'      ,...      
            'Fontsize'      ,16         ,...
            'FontWeight'    ,'Bold',...
            'FontName','Times',...
            'XMinorTick'    ,'on'       ,...
            'YMinorTick'    ,'on'       ,...
            'XColor'      , [0 0 0], ...
            'YColor'      , [0 0 0], ...
            'LineWidth'     ,0.5,...
            'Xgrid'         ,'off',...
            'Ygrid'         ,'off',... 
            'Box'           ,'on'...
        );
    set(gcf, 'color', 'white');


MeasurementNr = AVERAGES;
NumberOfMics = length(INPUT_CHAN)-1;
ReferenceMic = INPUT_CHAN==REFERENCE_MIC;
AvgF = zeros(length(Measurement.Nr1.freq),NumberOfMics+1);
for k=1:MeasurementNr
    MeasurementInstance = ['Nr',num2str(k)];
    for i=1:NumberOfMics
    H(:,i) = Measurement.(MeasurementInstance).Y(:,i)./Measurement.(MeasurementInstance).Y(:,ReferenceMic);
    end
    plot(AX3,Measurement.Nr1.freq(:,1),abs(H))
    plot(AX4,Measurement.Nr1.freq(:,1),angle(H))
    AvgF = AvgF + Measurement.(MeasurementInstance).Y./MeasurementNr ;
end
f_cal = Measurement.Nr1.freq(:,1);
for i=1:NumberOfMics
C(:,i) = AvgF(:,i) ./ AvgF(:,ReferenceMic);
end

for i=1:NumberOfMics
[PolReal(i,:),NonVal,MU_Real(i,:)] = polyfit(f_cal,real(C(:,i)),15);
[PolImag(i,:),NonVal,MU_Imag(i,:)] = polyfit(f_cal,imag(C(:,i)),15);
end

for i=1:NumberOfMics
    C_fitted(:,i) = polyval(PolReal(i,:),f_cal,[],MU_Real(i,:)) + 1i.*polyval(PolImag(i,:),f_cal,[],MU_Imag(i,:));
end
plot(AX3,f_cal,abs(C_fitted),'-.')
plot(AX4,f_cal,angle(C_fitted),'-.')
save(Data_File,'Measurement', 'f_cal', 'C','PolReal','PolImag','MU_Real','MU_Imag')

