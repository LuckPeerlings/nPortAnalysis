 function [] = VXI_TwoPort_Measurement()
%
% This function is used to control the VXI system to perform acoustic
% measurements using the low-level library
% 
close all
clc


% A. Defines  {{{1
%**************************************************************************
% A1. Defines that can be changed by the user ...

%A1.1 Important information is interactively asked

SAVEDIR = uigetdir; 
AMBIENT_TEMPERATURE_START = input('What is the ambient start temperature [Celsius]?  ');
AMBIENT_PRESSURE =          input('What is the ambient pressure  [hPa]?              ');
AMBIENT_HUMIDITY =          input('What is the relative humidity, in fraction  [-]?  ');
U =                         input('What is the flowspeed  [m/s]?                     ');
%        SAVEDIR = 'C:\Luck\TEST';
%        AMBIENT_TEMPERATURE_START = 20;
%        AMBIENT_PRESSURE =          1013.25;
%        AMBIENT_HUMIDITY =          0;
%        U =  0;

%A1.2 General measurement info

% See whether the user wants to use pre-saved data
reply =                    input('Do you want to use a configuration file? Y/N:     ', 's');
% reply = 'Y';
if strcmp(reply,'Y')
    [FILE,DIR]=uigetfile('.mat','Select configuration file');
    load([DIR,FILE]);
% load('H:\Luck\TEST\test.mat')
elseif strcmp(reply,'N')
    
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
else
    error('Invalid choice')
end
%A1.8 Check the inputs
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

%A2. Properties of the module.

CLOCK_FREQ = 51200;     % Module clock frequency
SPAN = 51200/(2.56*2^2);% Module span (used by inputs & for trigger timing)
% Allowed numbers for the span are: (see manual
% for more information)
% SPAN = CLOCKFREQ / (2.56*2^n), with n = 1,2,3,4.... OR
% SPAN = CLOCKFREQ / (2.56*5*2^n), with n = 1,2,3,4.... OR


BLOCKSIZE = 8192;     % Input Blocksize
BLOCKSIZE_S = 8192;   % Source Blocksize (also used in THIS demo to
%    determine the buffer size and length of data)
%    - Keep this an even number since buffersize is
%      set to BLOCKSIZE_S/2
%    - Keep this greater than 8192 for this demo
%      (this keeps the A & B buffers at least 4k in size)


TRIG_DELAY = 0;         % Input Trigger Delay in Input samples
% Pretriger: Minus values  Post Trigger: Plus Values

SPAN_S = 51200/2.56;    % Source span %Usually set to the highest value possible to have a smooth excitation signal



%A2.2 Check if the settings are corrcet
%Span Correct? See the manual for the allowed spans
if ~or((mod(log2(CLOCK_FREQ/(2.56*SPAN)),1) == 0),(mod(log2(CLOCK_FREQ*5/(2.56*SPAN)),1) == 0))
    warning('Input span is not allowed, device will set it to the nearest allowed value and the subsequent check on the blocksize may be invalid')
    disp('Press Enter to Continue')
    pause
    
end

%Does the blocksize contain an integer number of oscillations?
for ii = 1:length(MEASUREMENTFREQUENCY)
    if mod(BLOCKSIZE*MEASUREMENTFREQUENCY(ii)/(SPAN*2.56),1) ~= 0
        warning('The excitation Frequency %f, and possibly other frequencies do not create an integer number of oscillations in the measured blocksize ',MEASUREMENTFREQUENCY(ii))
        disp('Press Enter to Continue')
        pause
        break
    end
end

%**************************************************************************
%******************************************************************************

% A3. Other Defines

SOURCE_MODE = 'SOURCE_MODE_SINE';	% See VT1432 Help for possible modes
% Mapping arrays... use these to make user friendly disp statements.
% You can get the mappings from the VT1432A Help documentation (look
% under the help for the associated function)
SOURCE_MODES = {'SOURCE_MODE_SINE' 'SOURCE_MODE_BSINE' 'SOURCE_MODE_RAND' ...
    'SOURCE_MODE_BRAND' 'SOURCE_MODE_RANDZ' 'SOURCE_MODE_BRANDZ' ...
    'SOURCE_MODE_ARB' 'SOURCE_MODE_BARB'};
%}}}
% 1. Initialize Software connection to Hardware and Setup Channel Groups {{{1
StartDateTime = datestr(clock); %Save the starttime of the measurement
disp('1. Initializing Hardware...');
[status,session] = vt1432('init',ADDRESS,1,1);
vt1432_check_status(session,status);

% Find out about channels present and create appropriate groups of channels

[status,TotChans,NumInpChans,NumSrcChans,NumTachChans] = ...
    vt1432('getNumChans',session);
vt1432_check_status(session, status);

TotChans=double(TotChans);
NumInpChans=double(NumInpChans);
NumSrcChans=double(NumSrcChans);
NumTachChans=double(NumTachChans);

inp_chs_avail = 1:NumInpChans;						% Create list of input channel ids
Source_Chs_avail = 4097:4097+(NumSrcChans-1);	% Create list of source channel ids
tot_chs_avail = [inp_chs_avail Source_Chs_avail];	% Create list of all channel ids

Input_Chs = [MIC_CHAN REF_CHAN TEMP_CHAN ACCEL_CHAN CTA_CHAN OTHER_CHAN]; %Create list of the input channels
Sensitivity_Ref_Chan  = ones(size(REF_CHAN));
Sensitivity_Input_Chs = [SENSITIVITY_MIC_CHAN Sensitivity_Ref_Chan SENSITIVITY_TEMP_CHAN SENSITIVITY_ACCEL_CHAN SENSITIVITY_CTA_CHAN SENSITIVITY_OTHER_CHAN];
if length(Input_Chs) ~= length(Sensitivity_Input_Chs)
    error('The number of the channels and the number of sensitivities do not match')
end

Source_Chs =  OUTPUT_CHAN; %Create list of the source channels

% Check to make sure that all channels in Source_Chs and Input_Chs are
% actually present in the system

for i=1:length(Source_Chs)
    if isempty(find(Source_Chs_avail == Source_Chs(i),1))
        disp('Error: Source_Chs contains a channel id which is not present')
    end
end
disp(['     Source Channels On: ' num2str(Source_Chs)]);

for i=1:length(Input_Chs)
    if isempty(find(inp_chs_avail == Input_Chs(i),1))
        disp('Error: Input_Chs contains a channel id which is not present');
        disp('       Set Input_Chs = [] if there are no input channels');
    end
end
disp(['     Input Channels On: ' num2str(Input_Chs)]);

%create a group of all physical channels and then turn them off in case some
%previous activity left them on
[status,global_gid] = vt1432('createChannelGroup',session, ...
    length(tot_chs_avail),tot_chs_avail);
vt1432_check_status(session, status);
[status] = vt1432('setActive',session,global_gid,'CHANNEL_OFF');
vt1432_check_status(session, status);

%Create seperate Input
for ii = 1:length(Input_Chs)
    [status,gid_Single(ii)] = vt1432('createChannelGroup',session,1,Input_Chs(ii));
    vt1432_check_status(session,status);
end

%create seperate output
for ii = 1:length(Source_Chs)
    [status,gid_Single_Source(ii)] = vt1432('createChannelGroup',session,1,Source_Chs(ii));
    vt1432_check_status(session,status);
end



%create the seperate input-groups
[status,gid_Mic] = vt1432('createChannelGroup',session,length([MIC_CHAN]),[MIC_CHAN]);
vt1432_check_status(session,status);

[status,gid_Ref] = vt1432('createChannelGroup',session,length(REF_CHAN),REF_CHAN);
vt1432_check_status(session,status);

if ~isempty(MIC_CHAN_CALWAVE)
    Mic_Tot = [];
    for ss = 1:length(REF_CHAN)
        %Create the calibration signal groups for each port
        [status,gid_Mic_CalWave{ss}] = vt1432('createChannelGroup',session,length([MIC_CHAN_CALWAVE.(['Port',num2str(ss)]) REF_CHAN(ss)]),[MIC_CHAN_CALWAVE.(['Port',num2str(ss)]) REF_CHAN(ss)]);
        [MIC_CHAN_CALWAVE.(['Port',num2str(ss)]) REF_CHAN(ss)]
        vt1432_check_status(session,status);
        Mic_Tot = [Mic_Tot, MIC_CHAN_CALWAVE.(['Port',num2str(ss)])];
    end
    [status,gid_Mic_CalWave_Tot] = vt1432('createChannelGroup',session,length([Mic_Tot,REF_CHAN]),[Mic_Tot,REF_CHAN]);
    vt1432_check_status(session,status);
end
if ~isempty(TEMP_CHAN)
    [status,gid_Temp] = vt1432('createChannelGroup',session,length(TEMP_CHAN),TEMP_CHAN);
    vt1432_check_status(session,status);
end
if ~isempty(ACCEL_CHAN)
    [status,gid_ACCEL] = vt1432('createChannelGroup',session,length(ACCEL_CHAN),ACCEL_CHAN);
    vt1432_check_status(session,status);
end
if ~isempty(CTA_CHAN)
    [status,gid_CTA] = vt1432('createChannelGroup',session,length(CTA_CHAN),CTA_CHAN);
    vt1432_check_status(session,status);
end
if ~isempty(OTHER_CHAN)
    [status,gid_OTHER] = vt1432('createChannelGroup',session,length(OTHER_CHAN),OTHER_CHAN);
    vt1432_check_status(session,status);
end

%Create a global group of the inputs
[status,gid_i] = vt1432('createChannelGroup',session,length(Input_Chs),Input_Chs);
vt1432_check_status(session, status);

%Create a global group of the sources that will be active
[status,gid_s] = vt1432('createChannelGroup',session,length(Source_Chs),Source_Chs);
vt1432_check_status(session, status);



%Create a global group of inputs and sources
[status,gid_g] = vt1432('createChannelGroup',session,...
    length([Input_Chs Source_Chs]),[Input_Chs Source_Chs]);
vt1432_check_status(session, status);
%}}}
% 2. Setup "Global" parameters {{{1
disp('2. Setup "Global" parameters...');

[status] = vt1432('setClockFreq',session, gid_g, CLOCK_FREQ);
vt1432_check_status(session, status);
[status, gotClockFreq] = vt1432('getClockFreq', session, gid_g);
vt1432_check_status(session, status);
disp(['     Clock Freq: ' num2str(gotClockFreq)]);

% Source triggering of inputs is more accurate with DATA_SIZE_32
[status] = vt1432('setDataSize', session, gid_g, 'DATA_SIZE_32');
vt1432_check_status(session, status);

% Note: Span should be set even on VT1434 modules so that the source trigger
%       can line up with decimated input data (when other input modules are present).
[status] = vt1432('setSpan',session,gid_g,SPAN);
vt1432_check_status(session, status);
[status, Meas_Span] = vt1432('getSpan',session,gid_g);
vt1432_check_status(session, status);
disp(['     Span: ' num2str(Meas_Span)]);
%}}}
% 3. Setup Input Channel parameters (when input modules are present) {{{1

disp('3. Setup Input Channel parameters...');

[status] = vt1432('setActive',session,gid_i,'CHANNEL_ON');
vt1432_check_status(session, status);

[status] = vt1432('setBlocksize',session,gid_i,BLOCKSIZE);
vt1432_check_status(session, status);

[status, gotBlocksize] = vt1432('getBlocksize',session,gid_i);
vt1432_check_status(session, status);
disp(['     Input Blocksize: ' num2str(double(gotBlocksize))]);

% Set the range for each of the channels
Range_Chs = [MIC_RANGE REF_RANGE TEMP_RANGE ACCEL_RANGE CTA_RANGE OTHER_RANGE];
if length(Input_Chs) ~= length(Range_Chs)
    error('Length of the input channels is not equal to the length of the range vector')
end
for ii = 1:length(Input_Chs)
    [status] = vt1432('setRange',session,Input_Chs(ii),Range_Chs(ii));
    vt1432_check_status(session, status);
    [status, gotRange] = vt1432('getRange',session,Input_Chs(ii));
    vt1432_check_status(session, status);
    disp(['     Input Range of channel ',num2str(Input_Chs(ii)),': ', num2str(gotRange)]);
end
% Turn input AC coupling on to help remove any DC
[status] = vt1432('setCoupling', session, gid_i, 'COUPLING_AC');
vt1432_check_status(session, status);
% Turn input DC coupling on for the thermocouple signal, the CTA signal
% and the OTHER signals
if ~isempty(TEMP_CHAN)
    [status] = vt1432('setCoupling', session, gid_Temp, 'COUPLING_DC');
    vt1432_check_status(session, status);
end
if ~isempty(CTA_CHAN)
    [status] = vt1432('setCoupling', session, gid_CTA, 'COUPLING_DC');
    vt1432_check_status(session, status);
end
if ~isempty(OTHER_CHAN)
    [status] = vt1432('setCoupling', session, gid_OTHER, 'COUPLING_DC');
    vt1432_check_status(session, status);
end

%Set the ICP mode
if strcmp(ICP_MICROPHONE,'True')
    [status] = vt1432('setInputMode', session, gid_Mic, 'INPUT_MODE_ICP');
    vt1432_check_status(session, status);
    % After the trigger, there are some transient signals which need some time
    % to settle.
    pause(15)
end
% Input span already set above in global parms (so just get it for
% info purposes)
[status,SetSpan] = vt1432('getSpan',session,gid_i);
vt1432_check_status(session, status);
disp(['     Input span = ' num2str(SetSpan)]);

[status] = vt1432('setDataMode',session,gid_i,'BLOCK_MODE');
vt1432_check_status(session, status);


[status,maxBlockSize] = vt1432('getBlocksizeCurrentMax',session,gid_i);
vt1432_check_status(session, status);
disp(['     Max block size input = ' num2str(maxBlockSize)]);
%}}}
% 4. Setup Source Channel parameters {{{1

disp('4. Setup Source Channel parameters...');
% 4.1 Get and display source board info

% Get some source info and display it... this can be useful when reporting
% bugs or checking to make sure that all source boards have the same firmware
% revs in ROM.  Note: Sources must be turned on and the source buffer properly
% initialized for this function to work.  Also srcGetRev should not be called
% during an Arb measurement since it uses the source transfer buffers.
[status] = vt1432('setActive',session,gid_s,'CHANNEL_ON');
vt1432_check_status(session, status);


for i=1:length(Source_Chs),
    [status, romid, romdate, bdid, bddate] = vt1432('srcGetRev', session, ...
        Source_Chs(i));
    vt1432_check_status(session, status);
    disp(['     SOURCE CHAN: ' num2str(Source_Chs(i)) ...
        ' romid=' num2str(dec2hex(double(romid(1)))) ...
        ' romdate=' num2str(dec2hex(double(romdate(1)))) ...
        ' bdid=' num2str(dec2hex(double(bdid(1)))) ...
        ' bddate=' bddate]);
end

%4.2 Generic source parameters

[status] = vt1432('setSourceBlocksize',session, gid_s,BLOCKSIZE_S);
vt1432_check_status(session, status);
[status, gotSourceBlocksize] = vt1432('getSourceBlocksize', session, gid_s);
vt1432_check_status(session, status);
disp(['    Source Blocksize: ' num2str(double(gotSourceBlocksize))]);

[status] = vt1432('setRampRate',session,gid_s,RAMP_TIME_S);
vt1432_check_status(session, status);
[status, gotRampRate] = vt1432('getRampRate', session, gid_s);
vt1432_check_status(session, status);
disp(['     Source Ramp Rate: ' num2str(gotRampRate)]);

% Source uses a log DAC so pick the closest Range above AMP_S and
% then use AmpScale to scale the output down to the desired amplitude.
% (Output voltage = Range * AmpScale)
[status] = vt1432('setRange',session,gid_s,AMP_S);
vt1432_check_status(session, status);
[status, Source_Range] = vt1432('getRange',session,gid_s);
vt1432_check_status(session, status);
disp(['    Source Range: ' num2str(Source_Range)]);

[status] = vt1432('setAmpScale',session, gid_s, AMP_SCALE);
vt1432_check_status(session, status);
[status, gotAmpScale] = vt1432('getAmpScale',session,gid_s);
vt1432_check_status(session, status);
disp(['    Source Amp Scale: ' num2str(gotAmpScale)]);

[status] = vt1432('setSourceSpan', session, gid_s, SPAN_S);
vt1432_check_status(session, status);
[status, gotSourceSpan] = vt1432('getSourceSpan', session, gid_s);
vt1432_check_status(session, status);
disp(['    Source Span: ' num2str(gotSourceSpan) ' Hz']);

% 4.3 Source mode specific parameters
[status] = vt1432('setSourceMode',session,gid_s,SOURCE_MODE);
vt1432_check_status(session, status);
[status, gotSourceMode] = vt1432('getSourceMode', session, gid_s);
vt1432_check_status(session, status);
disp(['    Source Mode: ' SOURCE_MODES{double(gotSourceMode)-169}]);

[status] = vt1432('setSourceSpan', session, gid_s, SPAN_S);
vt1432_check_status(session, status);
[status, gotSourceSpan] = vt1432('getSourceSpan', session, gid_s);
vt1432_check_status(session, status);
disp(['    Source Span: ' num2str(gotSourceSpan) ' Hz']);

%}}}
% 5. Setup Triggering {{{1
disp('5. Setup Triggering...');

% Let all channels auto arm themselves
[status] = vt1432('setArmMode',session,gid_g,'AUTO_ARM');
vt1432_check_status(session, status);

% Let all channels auto trigger themselves
[status] = vt1432('setAutoTrigger',session,gid_g,'AUTO_TRIGGER');
vt1432_check_status(session, status);

% No trigger channels on since auto triggering is being used
[status] = vt1432('setTriggerChannel',session,gid_g,'CHANNEL_OFF');
vt1432_check_status(session, status);
%}}}
% 7. Start The Measurement {{{1

% Do an Auto Zero just before starting the measurement... shouldn't need to
% do an autozero during a measurement even when changing ranges.

%Turn the sources off (these do not have to be zeroed). This has to be
%done, otherwise the source will be turned on if the input channel and the
%source are on the same module
[status] = vt1432('setActive',session,gid_s,'CHANNEL_OFF');
vt1432_check_status(session, status);
%Zero the input channels
[status] = vt1432('autoZero',session,gid_i);
vt1432_check_status(session, status);

[status] = vt1432('setEnable',session,gid_i,'ENABLE_TYPE_TIME','ENABLE_OFF');
vt1432_check_status(session, status);


for ll = 1:NRREPETITIONS
    %Randomize the frequency vector
    if strcmp(SCRAMBLEFREQUENCY,'True')
            f = ScrambleFrequency(MEASUREMENTFREQUENCY,REF_CHAN);
    else
        for ii = 1:size(REF_CHAN,2)
            f(:,ii) = MEASUREMENTFREQUENCY;
        end
    end
    %If the amplitudes are obtained from a file load that file. 
    if strcmp(EXTERNAL_AMP_SCALE,'True') && ~strcmp(CALWAVE,'True');
        ExtData = load(DATAFILE);
        SaveFileAmpScale = ExtData.Measurement.Source.Scale;
    end
    %First perform the wave calibration and save the amplitudes, loop
    %over all the frequencies
    % 7a. Calibration of the input wave {{{2
    if strcmp(CALWAVE,'True') && ll == 1
        disp('7a. Starting the calibration of strength of the ingoing wave...');
        for ff = 1:size(f,1)
            %Creating the plots to reduce time.
            if ff == 1
                h_fig_CalPlot = figure;
                for ss = 1:size(REF_CHAN,2)
                    h_axes_CalPlot(ss) = subplot(1,size(REF_CHAN,2),ss);
                    cla(h_axes_CalPlot(ss),'reset')
                    hold(h_axes_CalPlot(ss),'on')
                    title(h_axes_CalPlot(ss),['Reference Source',num2str(ss),': ',num2str(f(ff,ss),'%5.0f'),' [Hz]'])
                    xlabel(h_axes_CalPlot(ss),'Number of tries')
                    ylabel(h_axes_CalPlot(ss),'[Pa] ')
                end
            else
                for ss = 1:size(REF_CHAN,2)
                    cla(h_axes_CalPlot(ss),'reset')
                    hold(h_axes_CalPlot(ss),'on')
                    title(h_axes_CalPlot(ss),['Reference Source',num2str(ss),': ',num2str(f(ff,ss),'%5.0f'),' [Hz]'])
                    xlabel(h_axes_CalPlot(ss),'Number of tries')
                    ylabel(h_axes_CalPlot(ss),'[Pa] ')
                end
            end
            %Initializng the channels
                
            for ss = 1:size(REF_CHAN,2)

                % Only enable those channels used for the measurement
                [status] = vt1432('setEnable',session,gid_Mic_CalWave{ss},'ENABLE_TYPE_TIME','ENABLE_ON');
                vt1432_check_status(session, status);
                [status] = vt1432('setActive',session,gid_Single_Source(ss),'CHANNEL_ON');
                vt1432_check_status(session, status);
                [status] = vt1432('setDataMode',session,gid_Mic_CalWave{ss},'CONTINUOUS_MODE');
                vt1432_check_status(session, status);
                % Set delay for the ramptime of the source
                [status] = vt1432('setTriggerDelay',session,gid_Mic_CalWave{ss},(SETTLETIME_CAL+RAMP_TIME_S)*2.56*SPAN);
                vt1432_check_status(session, status);
                % Setting frequencies
                [status] = vt1432('setSineFreq',session, gid_Single_Source(ss),f(ff,ss));
                vt1432_check_status(session, status);
                [status, SetSineFreq] = vt1432('getSineFreq', session, gid_Single_Source(ss));
                vt1432_check_status(session, status);
                disp(['    Source ',num2str(ss),' Sine Freq: ' num2str(SetSineFreq) 'Hz']);
                %Initializing the counters
                nn = 0; stop = 0; bb = 0; ref{ss} = []; signal{ss} = []; Success(ss) = false; data =[];
                NrChannels(ss)= length(MIC_CHAN_CALWAVE.(['Port',num2str(ss)]))+1;
                % HARD CODED THE SPEED OF SOUND
                c0 = sqrt((AMBIENT_TEMPERATURE_START+273.15)*1.401*287.058);
                A_inv{ss} = ModalMatrix(MIC_POS_CALWAVE.(['Port',num2str(ss)]),f(ff,ss),U*VELOCITY_CORR_CAL.(['Port',num2str(ss)]),c0);
            
                for qq = 1:length(MIC_CHAN_CALWAVE.(['Port',num2str(ss)]))
                    Sensitivity{ss}(qq) = SENSITIVITY_MIC_CHAN(MIC_CHAN == MIC_CHAN_CALWAVE.(['Port',num2str(ss)])(qq));
                end
            end
            % Starting the aquisition of the signals to perform the calibration of the
            % amplitude wave.
%             for ss=1:size(REF_CHAN,2)
%                [status] = vt1432('initMeasure',session,gid_Single_Source(ss));
%                vt1432_check_status(session, status);
%             end
            [status] = vt1432('initMeasure',session,gid_Mic_CalWave_Tot);
            vt1432_check_status(session, status);
            DataDiscarded = false;
            while ~stop
                %Discard the first
                %acquired data to remove the transients that can be
                %present.
                if ~DataDiscarded
                    DiscardData(SETTLETIME_CAL,session,...
                        gid_Mic_CalWave_Tot,length([Mic_Tot,REF_CHAN]),...
                        SPAN,BLOCKSIZE);

                    DataDiscarded = true;
                end
                
                [status, blockAvailable] =  vt1432('blockAvailable',session,global_gid); %Check if there is data available
                vt1432_check_status(session, status);
                if blockAvailable
                    if bb == 0
                        data = [];                        
                    end
                    %Obtaining the data for each source and saving it
                    %to different blocks
                    [status,data_r,count] = vt1432('readFloat64Data',session,gid_Mic_CalWave_Tot,'TIME_DATA', ...
                        BLOCKSIZE*sum(NrChannels),'WAIT_FLAG');
                    vt1432_check_status(session, status);
                    data = [data; reshape(data_r,BLOCKSIZE,sum(NrChannels))];  
                    bb = bb + 1; 
                end
                               
                if bb == NR_CAL_BLOCKS
                    %Reshape the data for each block
                    TotChan = 0;
                    for ss = 1:size(REF_CHAN,2)
                       signal{ss} = data(:,TotChan+1:TotChan+NrChannels(ss)-1);                        
                       TotChan = TotChan + NrChannels(ss)-1;                            
                       ref{ss} = data(:,end-size(REF_CHAN,2)+ss);
                    end
                    %Performing the wave decomposition for each source and
                    %scale the amplitude accordingly
                    bb = 0;
                    for ss = 1:size(REF_CHAN,2)
                        %Only do the calculation for those sources where
                        %the accuracy has not been reached yet.
                        if ~Success(ss)
                            %Scaling of the amplitude of the source
                            %Get the value of the amplitude scale
                            [status, OldAmpScale(ss)] = vt1432('getAmpScale',session,gid_Single_Source(ss));
                            vt1432_check_status(session, status);
                            %Calibrate the amplitudes and obtain the
                            %amplitude setting for the channel
                            [Success(ss), NewAmpScale(ss)] = CalibrateAmplitude(ref{ss},signal{ss},Sensitivity{ss},A_inv{ss}, OldAmpScale(ss),P_INGOING,CAL_TOL, h_axes_CalPlot(ss),nn);
                            % Update the amplitude scale
                            [status] = vt1432('setAmpScale',session, gid_Single_Source(ss), NewAmpScale(ss));
                            vt1432_check_status(session, status);
                            [status, NewAmpScale(ss)] = vt1432('getAmpScale',session,gid_Single_Source(ss));
                            vt1432_check_status(session, status);
                            disp(['Source ',num2str(ss),' Amp Scale : ' num2str(NewAmpScale(ss))]);
                            
                            %For the new iteration the first data has to be discarded
                            DataDiscarded = false;
                        end
                    end
                   nn = nn +1;
                   
                end
                if sum(Success == 1) == length(REF_CHAN)
                    disp('All the ports have reached their accuracy, breaking the calibration')
                    stop = 1;
                elseif nn == NR_TRYS
                    disp('    Maximum number of iterations reached for the calibration of the ingoing pressure wave')
                    stop = 1; 
                end                  
                
            end
            disp('    Closing the sources and input channels');
            for ss=1:size(REF_CHAN,2)
               [status] = vt1432('finishMeasure',session,gid_Single_Source(ss));
               vt1432_check_status(session, status);
            end
            [status] = vt1432('finishMeasure',session,gid_Mic_CalWave_Tot);
            vt1432_check_status(session, status);
            % Save the ampscale in the order of MEASUREMENTFREQUENCY
            for ss = 1:size(REF_CHAN,2)
                %Save the ampscale for the right frequency
                [status, SaveFileAmpScale(ss,f(ff,ss) == MEASUREMENTFREQUENCY)] = vt1432('getAmpScale',session,gid_Single_Source(ss));
                vt1432_check_status(session, status);
            end
        end
    end %}}}
    
    %Measurement of the data
    %Loop over the frequencies
    for ff = 1:size(f,1)
        %Make figure to show the uncertainty plot
        if ff == 1
            h_fig_UPlot = figure;
            for ss = 1:size(REF_CHAN,2)
                h_axes(ss) = subplot(1,size(REF_CHAN,2),ss);
                cla(h_axes(ss),'reset')
                hold(h_axes(ss),'on')
                title(h_axes(ss),['Reference Source',num2str(ss),': ',num2str(f(ff,ss),'%5.0f'),' [Hz]'])
                xlabel(h_axes(ss),'Number of measured blocks')
                ylabel(h_axes(ss),'Relative uncertainty [sqrt(sigma)/u]][%] ')
                set(h_axes,'YScale','log')
            end
        else
            for ss = 1:size(REF_CHAN,2)
                cla(h_axes(ss),'reset')
                hold(h_axes(ss),'on')
                title(h_axes(ss),['Reference Source',num2str(ss),': ',num2str(f(ff,ss),'%5.0f'),' [Hz]'])            
                xlabel(h_axes(ss),'Number of measured blocks')
                ylabel(h_axes(ss),'Relative uncertainty [sqrt(sigma)/u]][%] ')
                set(h_axes,'YScale','log')
            end
        end
        %Setting up the parameters for the data acquisition system
        %The VXI is set up to take continious data. If there is a
        %wait-time, blocks of data will be discarded. In this way, the
        %measurement only has to be initialized once.
        
        if ff == 1
            %Preallocating memory to speed up the measurement loop
            NrMeasBlocks = ceil(MEASUREMENTTIME*Meas_Span*2.56/BLOCKSIZE);
            Y = zeros(NrMeasBlocks*BLOCKSIZE,length(Input_Chs));
            NrExtraMeasBlocks = ceil(MAXEXTRATIME*Meas_Span*2.56/BLOCKSIZE);
            for ss = 1:size(REF_CHAN,2)
                MeanValue{ss} = zeros(length(gid_i),1);
                M2n{ss} = zeros(length(gid_i),1);
            end
            %Start the measurement, initializing sources
            [status] = vt1432('setDataMode',session,gid_i,'CONTINUOUS_MODE');
            vt1432_check_status(session, status);
            [status] = vt1432('setEnable',session,gid_i,'ENABLE_TYPE_TIME','ENABLE_ON');
            vt1432_check_status(session, status);
            %Set the triggerDelay equal to the ramp time of the source
            %and the waittime
            [status] = vt1432('setTriggerDelay',session,gid_i,(SETTLETIME+RAMP_TIME_S)*2.56*SPAN);
            vt1432_check_status(session, status);
        end
        %Setting the sine mode Frequency
        for ss = 1:length(Source_Chs)
            [status] = vt1432('setSineFreq',session, gid_Single_Source(ss),f(ff,ss));
            vt1432_check_status(session, status);
            [status, SetSineFreq] = vt1432('getSineFreq', session, gid_Single_Source(ss));
            vt1432_check_status(session, status);
            disp(['    Source ',num2str(ss),' Sine Freq: ' num2str(SetSineFreq) 'Hz']);
            %If the wave calibration has been performed, the source will be set to
            %the amplitude determined earlier.
            if strcmp(CALWAVE,'True') || strcmp(EXTERNAL_AMP_SCALE,'True')
                NewAmpScale = SaveFileAmpScale(ss,f(ff,ss) == MEASUREMENTFREQUENCY);
                [status] = vt1432('setAmpScale',session, gid_Single_Source(ss), NewAmpScale);
                vt1432_check_status(session, status);
                [status, NewAmpScale] = vt1432('getAmpScale',session,gid_Single_Source(ss));
                vt1432_check_status(session, status);
                disp(['Source ',num2str(ss),' Amp Scale : ' num2str(NewAmpScale)]);
            end
            [status] = vt1432('setActive',session,gid_Single_Source(ss),'CHANNEL_ON');
            vt1432_check_status(session, status);
        end
        
        % 7b. Performing the Measurement {{{2
        % In this measurement continuous data will be acquired
        
        
        [status] = vt1432('initMeasure',session,gid_s);
        vt1432_check_status(session, status);
        if strcmp(AUTORANGE,'True')
            disp('    AutoRanging...');
            for ii = 1:length(Input_Chs)
                [status] = vt1432('autoRange',session, gid_Single(ii),AUTORANGE_TIME);
                vt1432_check_status(session, status);
                [status,gotRange] = vt1432('getRange',session,gid_Single(ii));
                disp([' Range channel ', num2str(ii), ': ',num2str(gotRange)]);
                vt1432_check_status(session, status);
            end
        end
        [status] = vt1432('finishMeasure',session,gid_s);
        vt1432_check_status(session, status);
        
        disp('7b. Performing Measurement...');
        [status] = vt1432('initMeasure',session,global_gid);
        vt1432_check_status(session, status);

        %Enable in the input channels
        disp('    Acquiring data...');
        stop = 0;
        n = 0;
        m = 0;
        while ~stop
            [status, blockAvailable] =  vt1432('blockAvailable',session,gid_i); %Check if there is data available
            vt1432_check_status(session, status);
            if blockAvailable
                n = n+1;
                [status,data,count] = vt1432('readFloat64Data',session,gid_i,'TIME_DATA', ...
                    BLOCKSIZE*length(Input_Chs),'WAIT_FLAG');
                vt1432_check_status(session, status);
                data = reshape(data,BLOCKSIZE,length(Input_Chs));
                for ii = 1:length(Input_Chs)
                    Y((n-1)*BLOCKSIZE+1:n*BLOCKSIZE,ii) = data(:,ii);
                end
                
                for ss = 1:size(REF_CHAN,2)
                    [MeanValue{ss},M2n{ss},h_axes(ss),ExpectedUncertainty(:,ss)] = PlotUncertainty(h_axes(ss), REF_CHAN(ss), data, Input_Chs, MIC_CHAN,n,MeanValue{ss},M2n{ss});
                end
                for ss = 1:size(REF_CHAN,2)
                    g = sprintf('%06.2f ', MeanValue{ss}./SENSITIVITY_MIC_CHAN);
                    fprintf('Source %i: %s\n',ss,g)
                end
                if n > NrMeasBlocks
                    if strcmp(AUTODURATION,'True')
                        %Check if all the measurement channels have
                        %obtained the desired accuracy, otherwise
                        %measure longer.
                        if  m >= NrExtraMeasBlocks
                            fprintf('Maximum allowed extra time has been reached, stopping the measurement')
                            break
                            
                        elseif (sum(sum(ExpectedUncertainty > WANTEDACCURACY)) == 0)
                            stop = 1;
                            break
                        else
                            fprintf('Wanted Accuracy Not Obtained, One Extra Measurement Block Is Added')
                            m = m + 1;
                            continue
                        end
                    else
                        stop = 1;
                        break
                    end
                end
            end
        end
        %}}}
        
        disp('    Closing the sources and input channels');
        [status] = vt1432('finishMeasure',session,gid_i);
        vt1432_check_status(session, status);
        [status] = vt1432('finishMeasure',session,gid_s);
        vt1432_check_status(session, status);
        
        % 7c. Saving Time Data %{{{2
        
        disp('    Saving data...');
        savefile = ['Meas_',num2str(ff,'%03d'),'.mat'];
        savedir = [SAVEDIR,'\Rep_',num2str(ll)] ;
        if ~exist(savedir,'dir')
            mkdir(savedir)
        end
        save([savedir,'\',savefile],'Y');
        
        if strcmp(SCRAMBLEFREQUENCY, 'True')
            %Save the index of the frequency vectors to a temporary
            %file
            if ll ~= 1
                load([SAVEDIR,'\temp.mat'])
            end
            for ss = 1:size(REF_CHAN,2)
                [f_order,Index] = sort(f(:,ss));
                FreqIndex.(['Rep',num2str(ll,'%03i')])(:,ss) = Index;
            end
            save([SAVEDIR,'\temp.mat'],'FreqIndex')
        end
        
        if WAITTIME ~= 0;
            disp('    Waiting...')
            pause(WAITTIME)
        end
        
        
    end
end
%}}}
% 8a. Measurement of background signals {{{2
disp('8a. Measuring Background Signals');
if strcmp(MEASURE_BACKGROUND,'True');
    [status] = vt1432('setEnable',session,gid_i,'ENABLE_TYPE_TIME','ENABLE_OFF');
    vt1432_check_status(session, status);
    [status] = vt1432('setEnable',session,gid_i,'ENABLE_TYPE_TIME','ENABLE_ON');
    vt1432_check_status(session, status);
    [status] = vt1432('setDataMode',session,gid_i,'CONTINUOUS_MODE');
    vt1432_check_status(session, status);
    if strcmp(AUTORANGE_BACKGROUND,'True')
        disp('    AutoRanging...');
        for ii = 1:length(Input_Chs)
            [status] = vt1432('autoRange',session, gid_Single(ii),AUTORANGE_TIME);
            vt1432_check_status(session, status);
            [status,gotRange] = vt1432('getRange',session,gid_Single(ii));
            vt1432_check_status(session, status);
        end
    end
    disp('    Performing Background Measurement...');
    NrMeasBlocks = ceil(MEASUREMENTTIME_BACKGROUND*Meas_Span*2.56/BLOCKSIZE);
    Y = zeros(NrMeasBlocks*BLOCKSIZE,length(Input_Chs));
    % Start the measurement
    [status] = vt1432('initMeasure',session,gid_i);
    vt1432_check_status(session, status);
    n=0;
    stop = 0;
    while ~stop
        [status, blockAvailable] =  vt1432('blockAvailable',session,gid_i); %Check if there is data available
        vt1432_check_status(session, status);
        if blockAvailable
            [status,data,count] = vt1432('readFloat64Data',session,gid_i,'TIME_DATA', ...
                BLOCKSIZE*length(Input_Chs),'WAIT_FLAG');
            vt1432_check_status(session, status);
            data = reshape(data,BLOCKSIZE,length(Input_Chs));
            for ii = 1:length(Input_Chs)
                Y(n*BLOCKSIZE+1:(n+1)*BLOCKSIZE,ii) = data(:,ii);
            end
            n = n+1;
        end
        if n == NrMeasBlocks
            stop = 1;
        end
    end
    %Stop the measurement
    [status] = vt1432('finishMeasure',session,gid_s);
    vt1432_check_status(session, status);
    %Saving the background signals
    disp('    Saving data...');
    savefile = 'Background.mat';
    savedir = SAVEDIR;
    if ~exist(savedir,'dir')
        mkdir(savedir)
    end
    save([savedir,'\',savefile],'Y');
    
end
%}}}
% 9. Conclusion of the Measurement Loop
disp('9. Closing measurement system');
[status] = vt1432('setActive',session,global_gid,'CHANNEL_OFF');
vt1432_check_status(session, status);
[status] = vt1432('close',session);
vt1432_check_status(session, status);
%}}}
% 10. Creation of datafile to save {{{1
disp('10. Saving measurement data');
% Save the ambient data
Measurement.Date.Start = StartDateTime;
Measurement.Date.End = datestr(clock);
Measurement.RH = AMBIENT_HUMIDITY;
Measurement.t = AMBIENT_TEMPERATURE_START;
Measurement.p = AMBIENT_PRESSURE*100;
Measurement.U = U;
%Check if the measurement frequency is a row vector
if size(MEASUREMENTFREQUENCY,2) ~= 1
    MEASUREMENTFREQUENCY = MEASUREMENTFREQUENCY.';
end
Measurement.f = MEASUREMENTFREQUENCY;
Measurement.Repetitions = NRREPETITIONS;

Measurement.System.MeasTime = MEASUREMENTTIME;
Measurement.System.SampleFrequency = Meas_Span*2.56;
Measurement.System.BlockSize = BLOCKSIZE;

if strcmp(SCRAMBLEFREQUENCY, 'True')
    %Save the frequency indices to the measurementfile
    load([SAVEDIR,'\temp.mat'])
    Measurement.fIndex = FreqIndex;
else
    Measurement.fIndex = 1:length(MEASUREMENTFREQUENCY);
end

%If appropriate save the info of the wave calibration
if strcmp(CALWAVE, 'True')
    Measurement.WaveCal.Channel = MIC_CHAN_CALWAVE;
    Measurement.WaveCal.Velocity_Corr = VELOCITY_CORR_CAL;
    Measurement.WaveCal.Amplitude = P_INGOING;
    Measurement.WaveCal.Blocks = NR_CAL_BLOCKS;
end

%Saving the info on the channel-information
Measurement.Mic.Channel = MIC_CHAN;
Measurement.Mic.Pos = MIC_POSITION;
Measurement.Mic.Sensitivity = SENSITIVITY_MIC_CHAN;
Measurement.Ref.Channel = REF_CHAN;

if ~isempty(TEMP_CHAN)
    Measurement.Temp.Channel = TEMP_CHAN;
    Measurement.Temp.Sensitivity = SENSITIVITY_TEMP_CHAN;
    Measurement.Temp.ID = ID_TEMP_CHAN;
end
if ~isempty(ACCEL_CHAN)
    Measurement.Accel.Channel = ACCEL_CHAN;
    Measurement.Accel.Sensitivity = SENSITIVITY_ACCEL_CHAN;
    Measurement.Accel.ID = ID_ACCEL_CHAN;
    Measurement.Accel.Gain = GAIN_ACCEL_CHAN;
end
if ~isempty(CTA_CHAN)
    Measurement.CTA.Channel = CTA_CHAN;
    Measurement.CTA.Sensitivity = SENSITIVITY_CTA_CHAN;
    Measurement.CTA.ID = ID_CTA_CHAN;
end
if ~isempty(OTHER_CHAN)
    Measurement.Other.Channel = OTHER_CHAN;
    Measurement.Other.Sensitivity = SENSITIVITY_OTHER_CHAN;
    Measurement.Other.ID = ID_OTHER_CHAN;
end

% Saving information on the source
Measurement.Source.Range = AMP_S;
Measurement.Source.Scale = AMP_SCALE;
if strcmp(CALWAVE, 'True')
    Measurement.Source.Scale = SaveFileAmpScale;
end

%Saving information on the background measurement
if strcmp(MEASURE_BACKGROUND, 'True')
    Measurement.Background.Time = MEASUREMENTTIME_BACKGROUND;
end

savefile = 'MeasurementInfo.mat';
savedir = SAVEDIR;
if ~exist(savedir,'dir')
    mkdir(savedir)
end
save([savedir,'\',savefile],'Measurement');

% Example of the expected email structure:
% To be able to send a mail, get an account at mailgun.org and fill in the
% necessary values
email.APIkey  = '';
email.domain  = '';
email.from    = 'The seven MWL dwarves ';
email.to      = '';
email.subject = 'Measurements finished';
email.text    = 'Come down to the place without light, we are done!';

[status,answer] = mailgun(email);

end
%}}}

function DiscardData(Time,session,ID,NrChan,Span,BlockSize)
%Small Script to discard the acquired data.
BlocksToDiscard = ceil(Time*Span*2.56/BlockSize);
fprintf('Removing %2.1f sec of data. %i blocks will be discarded \n',Time,BlocksToDiscard)
fprintf('Blocks removed: ');

ii = 0;
while ii < BlocksToDiscard
    [status, blockAvailable] =  vt1432('blockAvailable',session,ID); %Check if there is data available
    vt1432_check_status(session, status);
    if blockAvailable
        [status,A,B] = vt1432('readFloat64Data',session,ID,'TIME_DATA',...
            BlockSize*NrChan,'WAIT_FLAG');
        vt1432_check_status(session, status);
        ii = ii + 1;
        fprintf('#')
    end
    pause(0.1)
end
fprintf('\n')
end
function f = ScrambleFrequency(MEASUREMENTFREQUENCY,REF_CHAN)
ii = 1;
f_list = MEASUREMENTFREQUENCY;
while isempty(f_list) == 0
    N = length(f_list);
    ii_rand = ceil(N*rand);%Generate random integer number
    f_rand(ii) = f_list(ii_rand);
    f_list(ii_rand) = [];
    ii = ii+1;
end
f(:,1) = f_rand;

%The other frequencies are compared against the
if size(REF_CHAN,2) == 2
    cc = 0; %Counter for the amount of tries to obtain an extra random vector
    ii = 1;
    f_list = MEASUREMENTFREQUENCY;
    while isempty(f_list) == 0
        N = length(f_list);
        ii_rand=ceil(N*rand);
        f_rand(ii) = f_list(ii_rand);
        if ( rem(f_rand(ii),f(ii)) == 0 || rem(f(ii),f_rand(ii)) == 0) ...
                || abs(f_rand(ii)/f(ii)-1) < 0.1
            disp('Frequency a multiple of each other or too close to each other')
            if N < 10
                disp('Not enough frequencies left, starting a new iteration')
                ii = 1;
                f_list = MEASUREMENTFREQUENCY;
                cc = cc + 1;
            end
        else
            f_list(ii_rand) = [];
            ii = ii+1;
        end
        if cc > 1000; error('Too many tries to find an appropriate frequency vector, check the settings'); end
    end
    f(:,2) = f_rand;
end
if size(REF_CHAN,2) > 2
    error('The script has not been coded to take care of more than two references')
end

end

function [MeanValue,M2n, h_axes, ExpectedUncertainty] = PlotUncertainty(h_axes, Ref_Chan, data, Input_Chs, Mic_Chan,n,MeanValuePrev,M2nPrev)
ColorList = colormap(lines(length(Mic_Chan)));
%Determining the SNR of the measurement and plotting this information
AnalyticSignal =  hilbert(data(:,Input_Chs == Ref_Chan)-mean(data(:,Input_Chs == Ref_Chan)));
AnalyticSignal =  AnalyticSignal./(mean(abs(AnalyticSignal)));
for ii = 1:length(Mic_Chan)
    PartialH(ii) = 2*conj(mean(AnalyticSignal.*data(:,ii)));
end
%Calculating the mean value and the unbiased variance
MeanValue = MeanValuePrev.*(n-1)./n + abs(PartialH)./n;
M2n = M2nPrev + (abs(PartialH) - MeanValuePrev).*(abs(PartialH)-MeanValue);

Sigma2 = M2n./(n);
ExpectedUncertainty = sqrt(Sigma2./n)./MeanValue;

for ii = 1:length(ExpectedUncertainty)
    plot(h_axes,n,ExpectedUncertainty(ii)*100,'o','markersize',10,'MarkerFaceColor',ColorList(ii,:),'MarkerEdgeColor',ColorList(ii,:));
end
drawnow
end

function [Success,NewAmpScale] = CalibrateAmplitude(ref,signal,Sensitivity,A_inv,OldAmpScale, P_ingoing,Cal_Tol,h_axes,nn)
%Performing the SynchronousDemodulation
%Creating the analytic signal
An_Signal = hilbert(ref);
%Creating a window
window = hann(length(ref));
for ii = 1:size(signal,2)
    %Multiply the analytic signal with the measured signal
    SigProduct = signal(:,ii)./An_Signal;
    %Taking the mean values (DC) value of the signal.
    %
    %The factor two is due to performing the demodulation and the value
    %has been scaled with the magnitude of the reference signal
    Z(ii) = 2*mean(SigProduct);
    %Applying the sensitivity to get to the value in Pascals
    Z(ii) = Z(ii)/Sensitivity(ii);
end
P_decomp = A_inv*Z.';
%Plotting the in and outgoing wave strengths
plot(h_axes,(nn+1),abs(P_decomp(2)),'x','markersize',10)
plot(h_axes,(nn+1),abs(P_decomp(1)),'o','markersize',10)

legend(h_axes,'Ingoing','Outgoing')
drawnow
P_measure = abs(P_decomp(2));

%Check if the ingoing wave is close to the measured wave,if so break
%the loop
Success = false;
if abs(P_ingoing/P_measure-1) < Cal_Tol
    disp('    Amplitude of the ingoing wave within the tolerance, breaking the calibration')
    NewAmpScale = OldAmpScale;
    Success = true;
    return
end

%Calculate the new amplitude scale, assuming a linear dependence
%between input and output
WantedAmpScale = OldAmpScale.*P_ingoing/P_measure;
NewAmpScale = WantedAmpScale;
if WantedAmpScale > 1
    NewAmpScale = 1;
end
%Check wether the change in signal amplitude is large and only allow a
%maximum change to avoid unstable behaviour
if abs(WantedAmpScale/OldAmpScale) < 0.01
    NewAmpScale = OldAmpScale*0.01;
    disp('    Warning: change in AmpScale is too large, NewAmpScale = OldAmpScale*0.75')
end
if abs(WantedAmpScale/OldAmpScale) > 100
    NewAmpScale = OldAmpScale*100;
    disp('    Warning: change in AmpScale is too large, NewAmpScale = OldAmpScale*1.75')
end
if NewAmpScale > 1
    NewAmpScale = 1;
    disp('    Warning: New amplitude scale is larger than one, amplitude scale is changed to 1')
end
if NewAmpScale == 1 && OldAmpScale == 1
    Success = true;
    disp(' Warning, the amplitude can not be increased any more, breaking the calibration')
    return
end
end

function Ainv = ModalMatrix(x,f,U,c0)

k_Up = 2*pi*f /c0 * 1/ (1-U/c0);
k_Down = 2*pi*f /c0 * 1/ (1+U/c0);
for ii = 1:length(x)
    A(ii,:) = [exp(-1i.*k_Down.*x(ii))  exp(1i.*k_Up.*x(ii))];
end
Ainv = inv(transp(A)*A)*transp(A);
end
