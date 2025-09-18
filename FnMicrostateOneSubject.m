function [MSResults,Preprocessed]=FnMicrostateOneSubject(EEGLabMat,Param)
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
%% -------------------- Loop over k results --------------------
for t=1:length(ClustRes.A_all)
    CurrentPrototypes=ClustRes.A_all{t};

    % Backfit
    [BackFittedLabels_IndividualMS,GMD]=MicroFit_GMD(EEGDataMat,CurrentPrototypes,Param.RespectPolarity);

    % Smooth
    [SmoothedBackFittedLabels_IndividualMS,RawLabels]=MicroSmooth_V2(EEGDataMat,CurrentPrototypes,Param.SmoothOpts.SmoothType,Param.SmoothOpts);

    % Statistics
    MicrostateStatSmooth=MicroStats(EEGDataMat,CurrentPrototypes,SmoothedBackFittedLabels_IndividualMS,Param.RespectPolarity,Preprocessed.fsample);
    MicrostateStat=MicroStats(EEGDataMat,CurrentPrototypes,BackFittedLabels_IndividualMS,Param.RespectPolarity,Preprocessed.fsample);

    %% -------------------- Store results (with index t) --------------------
    % GFP results (same for all k)
    MSResults.GFP.Avg=GFP.avg;
    MSResults.GFP.PeakIdxAll=PeakIdxAll;
    MSResults.GFP.NoiseIdx=PeakIdxAll(IsNoise);
    MSResults.GFP.PeakIdx=PeakIdx;

    % Store subject-specific microstate results per k
    MSResults.Prototypes{t}=CurrentPrototypes;
    MSResults.BackFittedLabels{t}=BackFittedLabels_IndividualMS;
    MSResults.SmoothedLabels{t}=SmoothedBackFittedLabels_IndividualMS;
    MSResults.StatsSmoothed{t}=MicrostateStatSmooth;
    MSResults.Stats{t}=MicrostateStat;

    % Store GMD
    MSResults.GMD{t}=GMD;

    % Clustering results
    MSResults.Modkmeans.Labels{t}=Labels;
    MSResults.Modkmeans.ClustRes{t}=ClustRes;
end
end
