function [] = DetermineComplexPressures()
DEFAULTCORRTIME = 10e-3;
CORR_THRESHOLD = 0.05;
MAXTIME = 5;
[DIRS] = uipickfiles();
[SAVEBASEDIR] = uigetdir('Choose the base directory to save the datafiles');

REPLY = input('Do you want to use determine the covariance matrix? Y/N ', 's');

if strcmp(REPLY,'Y')
    CalcCov = true;
elseif strcmp(REPLY,'N')
    CalcCov = false;
else
    error('Unknown input');
end

for dd = 1:length(DIRS)
    DIR = DIRS{dd};
    SAVEDIR = [SAVEBASEDIR,'\', fliplr(extractBefore(fliplr(DIR),'\'))]
    mkdir(SAVEDIR);
    disp(['Processing data in folder ', DIR])
    load([DIR,'\MeasurementInfo.mat']);
    
    if ~isfield(Measurement,'Repetitions')
        REPETITIONS = input('How many repetitions are there performed?');
    else
        REPETITIONS = Measurement.Repetitions;
    end
    
    for rr = 1:REPETITIONS
        
        DATADIR = [DIR,'\Rep_',num2str(rr)];
        
        %% Script to get the complex pressures from the measured time-data
        % All these variables have recopied from the workspace data.
        
        %Get the list of the data files in the
        list = dir(DATADIR);
        nn = 1;
        for ii = 1:length(list)
            if strncmpi(list(ii).name,'Meas',4)
                FILENAME = list(ii).name;
                FILENAME = regexprep(FILENAME, 'Meas_', '');
                FILENAME = regexprep(FILENAME, '.mat', '');
                MEASUREMENTNUMBER(nn) = str2double(FILENAME);
                nn = nn + 1;
            end
        end
        MEASUREMENTNUMBER = sort(MEASUREMENTNUMBER);
        MIC_CHAN = 1:length(Measurement.Mic.Channel);
        REF_CHAN = length(Measurement.Mic.Channel) + [1:size(Measurement.Ref.Channel,2)];
        
        if CalcCov
            if exist([DIR,'\Background.mat'],'file') == 0
                warning('The background signals have not been measured, reverting to a correlation time of %f s',DEFAULTCORRTIME)
                CorrTimeBackground = 0.01;
            else
                Background = load([DIR,'\Background.mat']);
                for ii = 1:length(Measurement.Mic.Channel)
                    CorrTimeBackground(ii) = CheckCorrelationTime(Background.Y(:,ii),Measurement.System.SampleFrequency,CORR_THRESHOLD,MAXTIME);
                    fprintf('Chan %i, Correlation time is %f [s] \n',ii,CorrTimeBackground(ii));
                end
            end
        end
        %% Getting the complex amplitudes wrt to the reference signal using the hilbert transform and filtering
        cc = 1;
        h = waitbar(cc/(size(REF_CHAN,2)*length(MEASUREMENTNUMBER)),'Wait');
        %Get all the frequencies of the references
        for ss = 1:size(REF_CHAN,2)
            for ff = 1:length(MEASUREMENTNUMBER)
                waitbar(cc/(size(REF_CHAN,2)*length(MEASUREMENTNUMBER)),h,'Wait');
                cc = cc +1;
                load([DATADIR,'\Meas_',num2str(round(MEASUREMENTNUMBER(ff)),'%03i'),'.mat'])
                
                ref = Y(:,REF_CHAN(ss)) - mean(Y(:,REF_CHAN(ss)));
                Ref_Freq(ss,ff) = GetRefFreq(ref,Measurement.System.SampleFrequency);
                MeasData.Ref_Freq(ff) = Ref_Freq(ss,ff);
                
                if CalcCov
                    %Get the windowsize for the excitation frequency and the plus
                    %and minus frequencies of the other
                    [windowSize,MovingAverageFilter,CorrTimeWindow] = GetFilterWindows(Ref_Freq(:,ff),Ref_Freq(ss,ff),Measurement.System.SampleFrequency,length(ref));
                    %Add the correlation time of the window to the set of correlation
                    %times.
                    CorrTime = [CorrTimeBackground,CorrTimeWindow];
                end
                
                %Normalizing the analytic signal, then the measured signals can be scaled directly using the sensitivity.
                %The absolute value of the reference is stored, so the transfer
                %function can be calculated back.
                MeasData.Ref(ff) = mean(abs(hilbert(ref)));
                an_signal = hilbert(ref)./MeasData.Ref(ff);
                
                for ii = 1:length(MIC_CHAN)
                    Proj = (Y(:,MIC_CHAN(ii))./Measurement.Mic.Sensitivity(ii))./an_signal;
                    
                    %Calculate the
                    if ~CalcCov
                        Z(ii) = mean(Proj);
                        continue
                    end
                    %Apply moving average filters
                    for tt = 1:length(windowSize)
                        Proj = filter(MovingAverageFilter{tt,:},1,Proj);
                        %Remove edge effects due to the moving average filter
                        
                    end
                    Proj(1:2*max(windowSize)) = [];
                    Proj(end-2*max(windowSize):end) = [];
                    
                    Z_SecondHarmonic(ii) = mean(Y(:,MIC_CHAN(ii))./Measurement.Mic.Sensitivity(ii)./an_signal.^2);
                    %Resample the projection according to the correlation time
                    %to have uncorrelated samples in the time series and store
                    %it in a 2D Vector to be able to calculate the
                    %cross-covariances between the other signals.
                    if ii == 1
                        clear Proj_Matrix
                    end
                    Proj_Matrix(ii,:) = Proj(1:ceil(max(CorrTime)*Measurement.System.SampleFrequency):end);
                    Z(ii) = mean(Proj_Matrix(ii,:));
                    if ii == 1
                        if length(Proj) < 100
                            warning('The length of the resampled projection is smaller than 100 samples, maybe problem with statistical data')
                        end
                    end
                end
                %Calculate the covariance
                if CalcCov
                    [CoVar, CompCoVar,CoVar_uu,CoVar_vv,CoVar_uv] = DetermineCoVariance(Proj_Matrix);
                end
                
                %Saving the variances
                MeasData.f(ff) = Measurement.f(Measurement.fIndex.(['Rep',num2str(rr,'%03i')])(ff,ss));
                MeasData.Z(:,ff) = Z;
                %Save also the other data when the covariance is calculated
                if CalcCov
                    MeasData.Z_SecondHarmonic(:,ff) = Z_SecondHarmonic;
                    MeasData.CoVar(ff,:,:) = CoVar; %Complex Covariance
                    MeasData.CompCoVar(ff,:,:) = CompCoVar; %Complex complementary Covariance
                    for ii = 1:length(MIC_CHAN)
                        MeasData.CoVar_NoCorr(ii,ff,1) = CoVar_uu(ii,ii);
                        MeasData.CoVar_NoCorr(ii,ff,2) = CoVar_uv(ii,ii);
                        MeasData.CoVar_NoCorr(ii,ff,3) = CoVar_uv(ii,ii);
                        MeasData.CoVar_NoCorr(ii,ff,4) = CoVar_vv(ii,ii);
                    end
                end
                if isfield(Measurement,'Accel')
                    Accel_Chan = length(Measurement.Mic.Channel) + length(Measurement.Ref.Channel) + [1:length(Measurement.Accel.Channel)];
                    for ii = 1:length(Accel_Chan)
                        
                        Proj = (Y(:,Accel_Chan(ii))./Measurement.Accel.Sensitivity(ii))./an_signal;
                        A(ii) = mean(Proj);
                    end
                    MeasData.Acceleration(:,ff) =  A;
                end
                
                %Saving the temperature data
                if isfield(Measurement,'Temp')
                    Temp_Chan = length(Measurement.Mic.Channel) + length(Measurement.Ref.Channel) + [1:length(Measurement.Temp.Channel)];
                    for ii=1:length(Measurement.Temp.Channel)
                        MeasData.T(ii,ff) = mean(Y(:,Temp_Chan(ii))) / Measurement.Temp.Sensitivity(ii) ;
                    end
                end
                if isfield(Measurement,'Other')
                    for ii=1:length(Measurement.Other.Channel)
                        MeasData.Other(ii,ff) = mean(Y(:,Measurement.Other.Channel(ii))) / Measurement.Other.Sensitivity(ii) ;
                    end
                end
                
            end
            if REPETITIONS == 1
                save([SAVEDIR,'\ComplexPressures_Source',num2str(ss),'.mat'],'MeasData');
                if ss == 1
                    if exist([SAVEDIR,'\MeasurementInfo.mat'],'file') == 2
                        warning('MeasurementInfo is not copied, it already exists in the folder')
                    else
                        copyfile([DIR,'\MeasurementInfo.mat'],SAVEDIR);
                    end
                end
            else
                if rr == 1 && ss==1
                    if exist([SAVEDIR,'\MeasurementInfo.mat'],'file') == 2
                        warning('MeasurementInfo is not copied, it already exists in the folder')
                    else
                        copyfile([DIR,'\MeasurementInfo.mat'],SAVEDIR);
                    end
                end
                save([SAVEDIR,'\Rep_',num2str(rr),'ComplexPressures_Source',num2str(ss),'.mat'],'MeasData');
                
            end
        end
        clear MeasData MEASUREMENTFREQUENCY
    end
    close(h)
end
end

function [CorrTime] = CheckCorrelationTime(Signal,Fs,Threshold,MaxTime)
%Determine the autocorrelation function and determine the position where it
%first crosses the determined treshold.
r = xcorr(Signal,Signal,MaxTime*Fs,'coeff');
ix = find(abs(r) > Threshold, 1, 'first');
CorrTime = (MaxTime*Fs + 1 - ix)/Fs;
end

function [Freq] = GetRefFreq(Signal,Fs)
[Spec,FreqVec] = cpsd(Signal,Signal,[],[],[],Fs);
[Val,I] = max(Spec);
Freq = FreqVec(I);
end

function [windowSize] = DetermineWindowSize(f,Fs)
if rem(Fs,f) ~= 0
    WindowFactor = rem(Fs,f)/f;
    [N,D] = rat(1/WindowFactor);
    %The fraction of a cycle that is missing to obtain a full period
    %Extending the windowSize to have an integer number of periods
    %in the averaging procedure
    warning('The filterfrequency is not a double of the excitationfrequency, increasing the window length by %i times',N)
    if N == 0
        N = 1;
    end
    windowSize = round(Fs/f*N);
else
    windowSize = round(Fs/f);
end
end

function [CoVar, CompCoVar,CoVar_uu,CoVar_vv,CoVar_uv] = DetermineCoVariance(Signal)
%Signal, first dimension the amount of channels, second dimension the
%temporal data.
CoVar_Single = zeros(size(Signal,1),size(Signal,1),size(Signal,2));
CompCoVar_Single = zeros(size(Signal,1),size(Signal,1),size(Signal,2));

MeanSignal = mean(Signal,2);
for ii = 1:size(Signal,2)
    Delta = Signal(:,ii)-MeanSignal;
    CoVar_Single(:,:,ii) = Delta*transpose(conj(Delta));
    CompCoVar_Single(:,:,ii) = Delta*transpose(Delta);
end

n = size(Signal,2);

% Unbiased variance of the sample mean
% Unbiased variance var_single = 1/(n-1) * sum x
% Variance of the mean var_mean = 1/n var_single
CoVar = 1/(n*(n-1))*sum(CoVar_Single,3);
CompCoVar = 1/(n*(n-1))*sum(CompCoVar_Single,3);

CoVar_uu = 1/2*real(CoVar + CompCoVar);
CoVar_vv = 1/2*real(CoVar - CompCoVar);
CoVar_uv = 1/2*imag(CompCoVar - CoVar);

end

function [windowSize,MovingAverageFilter,CorrTime] = GetFilterWindows(SourceFrequencies,Ref_Freq,Fs,SignalLength)
cc = 1;
for ii = 1:length(SourceFrequencies)
    if SourceFrequencies(ii) ~= Ref_Freq
        %Get the differences of the other frequencies
        FilterFreq(cc) = SourceFrequencies(ii) + Ref_Freq;
        FilterFreq(cc+1) = abs(SourceFrequencies(ii) - Ref_Freq);
        cc = cc +2;
    else
        %If the filter frequencies is the same as the excitation frequency,
        %the double frequency has to be filtered.
        FilterFreq(cc) = Ref_Freq*2;
        cc = cc + 1;
    end
end
for ii = 1:length(FilterFreq)
    windowSize(ii) = DetermineWindowSize(FilterFreq(ii),Fs);
    MovingAverageFilter{ii,:} = (1/windowSize(ii))*ones(1,windowSize(ii));
    CorrTime(ii) = windowSize(ii)/Fs;
    if windowSize(ii)/SignalLength > 0.01;
        warning('Window size, %i, is larger then 0.01 of the the complete signal length',windowSize(ii))
    end
end
end