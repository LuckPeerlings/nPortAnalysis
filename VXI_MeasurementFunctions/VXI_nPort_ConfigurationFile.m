% Script to create the configuration file for the two-port VXI
% measurement system.
close all
clear all
[FILE,DIR]= uiputfile('.mat','Select the destination file');

%The address of the VXI-system
    ADDRESS =  'VXI0::15,8::INSTR';
    % How long should the signals be aquired for each frequency?
    MEASUREMENTTIME = 20;
    % Which frequencies should be measured?
    MEASUREMENTFREQUENCY =  [2000:100:2700].';
    %Time to wait between measurements
    WAITTIME = 0.1;
    %The amount of repetitions one wants to have.
    NRREPETITIONS = 1;
    
    %Set if the measuremenent time has to be auto-adjusted such that the
    %wanted accuracy is reached
    AUTODURATION = 'False';
    WANTEDACCURACY =  5; % The wanted accuracy is defined as sigma/mu and is an estimate of the relative uncertainty
    MAXEXTRATIME = 5; %The maximum extra time (seconds) that will be added to the measurement
    
    
    %A1.3 Input Channels
    % list of input channels grouped by measurement type to turn on...
    %    - set this to [] if you do not have or want to turn on any input channels
    %    - this can be a sparse array (such as [1 3 5])
    %    - also add the logical address of any cards to the vt1432_init
    %      line later in this code.
    %    - Channel numbering should be in ascending order
    
    MIC_CHAN = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16];
    ICP_MICROPHONE = 'False'; %The breakout box only accepts blocks of 4 channels to be set to ICP at the same time
    %If there is multiple reference signals added, the program will excite
    %each of the reference at the same time, with frequencies that will not
    %overlap using the current blocksize.
    REF_CHAN = [17,18];
    TEMP_CHAN = [];
    ACCEL_CHAN = [];
    CTA_CHAN = [];
    OTHER_CHAN = [];
    
    % Set the range of each of the channel manually
    % If you set the switch auto range to true these settings will be
    % overridden
    MIC_RANGE = [1,1,1,1,...
                 1,1,1,1,...
                 1,1,1,1,...
                 1,1,1,1];
    REF_RANGE = [5,5];
    TEMP_RANGE = [];
    ACCEL_RANGE =[];
    CTA_RANGE =[];
    OTHER_RANGE =[];
    
    % Set the sensitivity of each of the channels. This information will be
    % stored in the measurement file, such that one can easily convert to
    % engineering units.
    % This information is also used during the adjustment of the ingoing
    % wave strength, if this option is chosen.
    SENSITIVITY_MIC_CHAN = [1,1,1,1,...
                            1,1,1,1,...
                            1,1,1,1,...
                            1,1,1,1]*31.6e-3;
      
    
    SENSITIVITY_REF_CHAN = [1,1];
    SENSITIVITY_TEMP_CHAN = [];
    SENSITIVITY_ACCEL_CHAN = [];
    SENSITIVITY_CTA_CHAN = [];
    SENSITIVITY_OTHER_CHAN = [];
    
    GAIN_ACCEL_CHAN = [];
    
    % Indentification number for the sensors used in the channels. This
    % information will be saved in the measurement file
    ID_MIC_CHAN = { '0011','0009','0001','0022',...
                    '0014','0016','0008','0006',...
                    '0020','0023','0002','0018',...
                    '0012','0017','0010','0021'...
                  };   
    ID_REF_CHAN = {'x','x'};
    ID_TEMP_CHAN = {};
    ID_ACCEL_CHAN = {};
    ID_CTA_CHAN = {};
    ID_OTHER_CHAN = {};
    
    %The microphone positions for each of the ports
    MIC_POSITION.Port1 = [  0, 195, 250, 400,...
                            455, 510, 565, 620,...
                            675, 730, 785, 840,...
                            895, 1100, 1155, 1350 ...
                          ]*1e-3;    
    
    %A1.4 Wave Calibration
    
    %Set if the wave calibration should be performed.
    CALWAVE = 'True';
    %Set if an external file should be used to set the strength of the
    %source. When using the external amp scale, the wave calibration is not
    %performed, but the strength of the source is obtained from the file.
    %This is helpful in the case when the system under study does not
    %change (a lot) based from a previous file. The wave calibration in
    %this case is skipped, significantlly reducing the measurement time.
    EXTERNAL_AMP_SCALE = 'False';
    DATAFILE = 'H:\Luck\RigidPlate2017_04_10\Sand_TwoMicrophones\MeasurementInfo.mat';
    %Parameters for the wave calibration, only one block of data is used to
    %determine the amplitude of the wave.
    
    %For each port, the information for the wave calibration should be
    %set. The measurement function will excite the corresponding source,
    %that is for Port1, the first element in the OUTPUT_CHAN vector is used
    %as source. For Port2, the second element is used and so on.
    
    %The information for each of the ports should be placed in a structure
    %with field name Port1, Port2 and so on.
    %Channels to use for the calibration, they have to be in ascending
    %order and the reference should have a higher channel number then the
    %microphone signals.
    MIC_CHAN_CALWAVE.Port1 = [1,2,3];
    MIC_POS_CALWAVE.Port1 = [250, 55, 0]*1e-3;
    %Velocity correction to determine the mean velocity in the measurement
    %port. This velocity correction is used to get the correct mean flow
    %velocity.
    VELOCITY_CORR_CAL.Port1 = 0.8978;
    
    %Information for the second port    
    MIC_CHAN_CALWAVE.Port2 = [14,15,16];
    MIC_POS_CALWAVE.Port2 = [0, 55, 250]*1e-3;
    VELOCITY_CORR_CAL.Port2 = -0.8978;
    
    %Maximum number of iterations to find the right ingoing wave strength
    NR_TRYS = 10;
    %Strength of the ingoing wave
    P_INGOING = 1; %[Pa]
    %Number of blocks to used to calibrate the wave. The more blocks, the
    %more time data is used for each iteration and the SNR increases. The
    %downside is that the calibration procedure takes longer time.
    NR_CAL_BLOCKS = 8;
    % How close should the measured value be to the wanted value, if the
    % relative difference is lower than CAL_TOL the calibration is finished
    CAL_TOL = 0.05;
    %Time to wait after adjusting the source amplitude, before acquiring
    %the data.
    SETTLETIME_CAL = 0.2;
    
    
    
    %A1.5 Source parameters
    % list of source channels to use for the excitation
    OUTPUT_CHAN = [1,2] + 4096;
    % Desired output level of the source is equal to AMP_S*AMP_SCALE
    % Range of the output channel
    AMP_S = 5; %[V]
    % Scale of the Amplitude
    AMP_SCALE = 0.05;
    RAMP_TIME_S = 0.1;      % Source ramp time in seconds
    %Time to wait after the ramp time of the source to start measuring
    SETTLETIME = 0.2;
    
    
    %A1.6 Background measurement
    % Should the background be measured
    MEASURE_BACKGROUND = 'False';
    AUTORANGE_BACKGROUND = 'False';
    MEASUREMENTTIME_BACKGROUND = 10;
    
    %A1.7 Other parameters
    %Should the frequencies be scrambled during the measurement?
    SCRAMBLEFREQUENCY = 'True';
    %Should the input channels be autoranged before each measurement
    AUTORANGE = 'True';
    %Time needed for the autoranging [s]
    AUTORANGE_TIME = 0.1;
   
%A1.8 Check the inputs
%     if length(OUTPUT_CHAN) ~= length(SPEAKER_POSITION)
%        error('Length of the speaker positions is not equal to the number of channels')
%     end
if strcmp(EXTERNAL_AMP_SCALE,'True') && strcmp(CALWAVE,'True')
    warning('Both the wave calibration and the amplitude obtained from a file are set to true')
end
if length(MEASUREMENTFREQUENCY) > 999
    error('Too many frequencies, > 999, problem with the naming of the savefiles')
end
if length(REF_CHAN) ~= length(OUTPUT_CHAN)
    error('The number of reference channels is not equal to the output channels')
end
if ~((length(MIC_CHAN) == length(SENSITIVITY_MIC_CHAN)) && (length(MIC_CHAN) == length(ID_MIC_CHAN)))
    error('Microphone inputs do not have the same length')
end
if ~((length(TEMP_CHAN) == length(SENSITIVITY_TEMP_CHAN)) && (length(TEMP_CHAN) == length(ID_TEMP_CHAN)))
    error('Temperature inputs do not have complete information')
end
if ~((length(ACCEL_CHAN) == length(SENSITIVITY_ACCEL_CHAN)) && (length(ACCEL_CHAN) == length(ID_ACCEL_CHAN)))
    error('Accelerometer inputs do not have complete information')
end
if ~((length(CTA_CHAN) == length(SENSITIVITY_CTA_CHAN)) && (length(CTA_CHAN) == length(ID_CTA_CHAN)))
    error('CTA inputs do not have complete information')
end
if ~((length(OTHER_CHAN) == length(SENSITIVITY_OTHER_CHAN)) && (length(OTHER_CHAN) == length(ID_OTHER_CHAN)))
    error('The extra inputs do not have complete information')
end

fn = fieldnames(MIC_POSITION);
nn = 0;
for ii = 1:length(fn)
    nn = nn + length(MIC_POSITION.(fn{ii}));
end
if nn ~= length(MIC_CHAN)
    error('The number of microphone positions does not match with the number of microphone channels')
end
for ii = 1:length(fn)
    if ~isfield(MIC_CHAN_CALWAVE,fn{ii});
        error('The fields in the MIC_POSITION do not correspond with the fields in MIC_CHAN_CALWAVE')
    end
    if length(MIC_CHAN_CALWAVE.(fn{ii})) > length(MIC_POSITION.(fn{ii}))
        error('The amount of microphones in the calibration is larger the avaliable microphone positions')
    end
end
clear fn nn

save([DIR,FILE],'-regexp', '^(?!(DIR|FILE)$).')