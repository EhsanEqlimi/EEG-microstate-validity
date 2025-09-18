function Param=FnParamMaker4Microstate()
% FnParamMaker_Microstates - Create a parameter structure for EEG microstate analysis
%
% OUTPUT:
%   Param - structure containing all settings for microstate pipeline
%
% Author: Ehsan Eqlimi, 30/11/2024, WAVES group, Ghent University, Belgium
% Email : Ehsan.Eqlimi@ugent.be, Ehsan.Eqlimi@outlook.com
% Updated: 2025 to support full pipeline (clustering + backfitting + smoothing + stats)

%% ----------------- Saving & Visualization -----------------
Param.DoSaveFolder=0;            % 0: do not auto-create subject folders, 1: create folders
Param.SaveMSResults=1;           % Save microstate clustering results
Param.SaveAllNonSortedMSFigs=1;  % Save all microstate maps before sorting
Param.PopoUpNonSortedMsFigs=1;   % Display non-sorted microstate figures
Param.PopoUpSortedMsFigs=1;      % Display sorted microstate figures
Param.SaveAllSortedMSFigs=1;     % Save all sorted microstate figures

%% ----------------- Preprocessing -----------------
Param.DoReRef='yes';             % Re-reference EEG to average reference
Param.DoDeTrend='yes';           % Remove linear trend from EEG
Param.DoDeMean='yes';            % Remove mean from EEG channels
Param.DoBandPassFilter='yes';    % Apply band-pass filter
Param.BandPassFrequencyBand=[1,45]; % Frequency range for band-pass (Hz)
Param.BandPassFilterType='firws';   % Filter type (finite impulse response, windowed-sinc)

%% ----------------- GFP (Global Field Power) settings -----------------
Param.GFPCalcMethod='amplitude'; % Method to compute GFP (amplitude or variance)
Param.DoNormalizeGFP=1;          % Normalize GFP maps before clustering
Param.MinPeakDist=10;            % Minimum distance (ms) between GFP peaks
Param.GFPthresh=1;               % Threshold for GFP peak selection (in std units)

%% ----------------- Clustering -----------------
Param.KRange=4;%:10;               % Range of microstate numbers to test
Param.Opts.reps=100;             % Number of clustering repetitions for stability
Param.Opts.max_iterations=1000;  % Maximum iterations per clustering run
Param.Opts.thresh=1e-6;          % Convergence threshold for clustering
Param.Opts.fitmeas='GEV';        % Goodness-of-fit measure: 'GEV' (Global Explained Variance) or 'CV'
Param.Opts.optimised=1;          % Whether to optimize cluster assignment

%% ----------------- Sorting -----------------
Param.RespectPolarity=0;         % 0=ignore polarity, 1=respect polarity of microstates

%% ----------------- Backfitting + Smoothing -----------------
Param.BackfitMethod='global';    % Backfitting method: 'global' (use group template) or 'subject-specific'

Param.SmoothOpts.b=3;            % Number of neighboring time points for smoothing
Param.SmoothOpts.lambda=5;       % Regularization parameter for smoothing
Param.SmoothOpts.max_iterations=1000; % Max iterations for label smoothing
Param.SmoothOpts.thresh=1e-6;    % Threshold for convergence in smoothing
Param.SmoothOpts.SmoothType='reject segments'; % Smoothing approach ('reject segments' or 'windowed')
Param.SmoothOpts.minTime=30;     % Minimum segment duration in ms for smoothing

%% ----------------- Visualization -----------------
Param.DoRdBuColorMap=1;          % Use red-blue colormap for topographical plots
end