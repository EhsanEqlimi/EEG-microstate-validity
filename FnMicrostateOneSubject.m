function MSResults=FnMicrostateOneSubject(EEGLabMat,Param)
% RunMicrostateOneSubject - Process EEG data for one subject
%
% INPUTS:
%   EEGLabMat      : struct loaded from a .mat file containing EEGLAB struct (ALLEEG or EEG)
%   TemplateFile   : .mat file with predefined microstate template (struct MSTemplate.Prototypes)
%   ChanLocFile    : channel location file (.ced or .locs)
%   ResultPath     : output folder for saving results
%   Param          : parameter struct from FnParamMaker4MS4()
%
% OUTPUT:
%   MSResults      : struct with individualized and predefined-template backfitting results
%
% Author: Ehsan Eqlimi, PhD @Ghent University, BE
% Email : ehsan.eqlimi@outlook.com
% Date  : 2025-09-18

%% -------------------- Extract EEGLAB struct --------------------
% Check if ALLEEG field exists and assign, otherwise search manually
if isfield(EEGLabMat,'ALLEEG')
    ALLEEG=EEGLabMat.ALLEEG;          % Directly assign ALLEEG if present
else
    ALLEEG=[];
    Candidate=EEGLabMat;             % Otherwise, take the main variable
    if isstruct(Candidate)&&isfield(Candidate,'chanlocs')&&isfield(Candidate,'data')
        ALLEEG=Candidate;            % Assign if it has EEG-like structure
    end  
end
if isempty(ALLEEG)
    error('No valid EEGLAB struct found in input matFile'); % Error if no EEG struct found
end
if iscell(ALLEEG),ALLEEG=ALLEEG{1};end  % Handle cell arrays
if length(ALLEEG)>1,ALLEEG=ALLEEG(1);end  % Use first element if multiple structs

%% -------------------- Convert to FieldTrip --------------------
EEGData=eeglab2fieldtrip(ALLEEG,'raw','none');  % Convert EEGLAB struct to FieldTrip format
if isfield(EEGData,'elec')&&isfield(EEGData.elec,'elecpos')
    EEGData.elec.chanpos=EEGData.elec.elecpos;    % Ensure electrode positions exist
end
EEGData.unit='mm';                               % Set unit for electrode positions

% Configure preprocessing: re-reference, demean, detrend
cfg=[];
cfg.reref=Param.DoReRef;                         % Re-reference if requested
cfg.refchannel='all';                             % Use all channels for referencing
cfg.demean=Param.DoDeMean;                        % Remove mean
cfg.detrend=Param.DoDeTrend;                      % Remove linear trend
Preprocessed=ft_preprocessing(cfg,EEGData);      % Apply preprocessing

%  Band-pass filtering
if Param.DoBandPassFilter
    cfg=[];
    cfg.bpfilter=Param.DoBandPassFilter;         % Enable band-pass filter
    cfg.bpfilttype=Param.BandPassFilterType;     % Filter type (e.g., FIR/IIR)
    cfg.bpfreq=Param.BandPassFrequencyBand;      % Frequency range
    Preprocessed=ft_preprocessing(cfg,Preprocessed); % Apply band-pass filter
end

%% -------------------- GFP + clustering --------------------
cfg=[];
cfg.method=Param.GFPCalcMethod;                  % Choose GFP calculation method
GFP=ft_globalmeanfield(cfg,Preprocessed);       % Compute Global Field Power (GFP)

% Detect GFP peaks
MinTFdistance=Param.MinPeakDist*Preprocessed.fsample/1000;  % Convert ms to samples
[~,PeakIdxAll]=findpeaks(GFP.avg,'MinPeakDistance',MinTFdistance); % Find peaks
IsNoise=GFP.avg(PeakIdxAll)>(mean(GFP.avg)+Param.GFPthresh*std(GFP.avg)); % Thresholding
PeakIdx=PeakIdxAll(~IsNoise);                   % Keep only non-noise peaks

% Extract EEG samples at GFP peaks
EEGDataMat=cell2mat(Preprocessed.trial);           % Concatenate trials
GFPData=EEGDataMat(:,PeakIdx);                     % Take columns at GFP peaks
if Param.DoNormalizeGFP
    GFPData=GFPData./mean(std(GFPData,0,2));   % Normalize GFP data across channels
end

% Apply modified k-means to identify microstate prototypes
[Prototypes,Labels,ClustRes]=modkmeans(GFPData,Param.KRange,Param.Opts);

%% -------------------- Backfitting to own microstates --------------------
% Backfit the EEG data using only the subject's own individualized microstates
% No grand-average template or external templates are used

% Backfit using the microstate prototypes obtained from this subject
[BackFittedLabels_IndividualMS,GMD]=MicroFit(EEGDataMat,Prototypes,Param.RespectPolarity);
%% -------------------- Smooth the backfitted labels -------------------
% Apply smoothing using the prototypes
SmoothedBackFittedLabels_IndividualMS=MicroSmooth(EEGDataMat,Prototypes,Param.SmoothOpts.SmoothType,Param.SmoothOpts);

%% -------------------- Statistics --------------------
% Compute microstate statistics
MicrostateStatSmooth=MicroStats(EEGDataMat,Prototypes,SmoothedBackFittedLabels_IndividualMS,Param.RespectPolarity,EEGDataMat.fsample);
MicrostateStat=MicroStats(EEGDataMat,Prototypes,BackFittedLabels_IndividualMS,Param.RespectPolarity,EEGDataMat.fsample);
%% -------------------- Return results --------------------
% Store GFP-related results
MSResults.GFP.Avg=GFP.avg;                         % Global Field Power time series
MSResults.GFP.PeakIdxAll=PeakIdxAll;               % All detected GFP peaks
MSResults.GFP.NoiseIdx=PeakIdxAll(IsNoise);        % Peaks considered as noise
MSResults.GFP.PeakIdx=PeakIdx;                     % Final peaks used for clustering

% Store subject-specific microstate results
MSResults.Prototypes=Prototypes;                  % Microstate topographies (maps)
MSResults.BackFittedLabels=BackFittedLabels_IndividualMS; % Labels assigned before smoothing
MSResults.SmoothedLabels=SmoothedBackFittedLabels_IndividualMS; % Labels after smoothing
MSResults.StatsSmoothed=MicrostateStatSmooth;    % Microstate statistics for smoothed labels
MSResults.Stats=MicrostateStat;                  % Microstate statistics for raw (unsmoothed) labels

% Store Global Map Dissimilarity (GMD) before smoothing
MSResults.GMD=GMD;                                % Measure of map dissimilarity

% Store k-means clustering results
MSResults.Modkmeans.Labels=Labels;               % Cluster assignments for GFP peaks
MSResults.Modkmeans.ClustRes=ClustRes;           % Detailed clustering results (e.g., centroids, GEV)

end
