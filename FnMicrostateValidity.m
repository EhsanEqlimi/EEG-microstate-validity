function MSResults = FnMicrostateValidity(Preprocessed, Param)
%% -------------------- Microstate Validity Indicator (MVI) --------------------
% Computes global and time-varying MVI for EEG data
% Inputs:
%   Preprocessed: FieldTrip-style EEG data (Preprocessed.trial, Preprocessed.fsample)
%   Param (optional): structure, e.g. Param.Bands = [4 7; 8 12; 13 30; 31 45];
% Output:
%   MSResults: struct containing global and time-varying MVI metrics

if nargin < 2
    Param = struct();
end
if ~isfield(Param,'Bands')
    Param.Bands = [4 7; 8 12; 13 30; 31 45]; % theta, alpha, beta, gamma
end

EEGDataMat = cell2mat(Preprocessed.trial);   % nChannels x nSamples
fs = Preprocessed.fsample;

%% --- Global covariance
Cx = cov(EEGDataMat');
lambda = sort(eig(Cx),'descend');
FE1   = lambda(1)/sum(lambda);
reff  = (sum(lambda))^2 / sum(lambda.^2);

%% --- Band-specific covariances
numBands = size(Param.Bands,1);
phi_b = zeros(1,numBands);
for b = 1:numBands
    band = Param.Bands(b,:);
    try
        cfg = [];
        cfg.bpfilter = 'yes';
        cfg.bpfreq = band;
        tmp = ft_preprocessing(cfg, Preprocessed);
        Cb  = cov(cell2mat(tmp.trial)');
    catch
        % fallback filtering
        nyq  = fs/2;
        W    = band/nyq;
        b_filt = fir1(256,W,'bandpass');
        data_f = filtfilt(b_filt,1,double(EEGDataMat'))';
        Cb = cov(data_f');
    end
    phi_b(b) = trace(Cb)/trace(Cx);
end
phi_max = max(phi_b);

%% --- Rényi Information Dimension proxy (global)
dR = estimateRID(EEGDataMat, fs);

%% --- Global MVI
MVI = (FE1/reff) * (1/(1+(dR-1))) * phi_max;
if MVI > 0.5
    ValidityDecision = 'One-hot microstate model valid';
else
    ValidityDecision = 'One-hot model likely invalid; use subspace-aware methods';
end

%% -------------------- Time-varying features --------------------
winMs   = 200; stepMs = 50;  % window size and step
winSamps  = round(winMs*fs/1000);
stepSamps = round(stepMs*fs/1000);
nWins = floor((size(EEGDataMat,2)-winSamps)/stepSamps)+1;

FE1_tv   = zeros(1,nWins);
reff_tv  = zeros(1,nWins);
dR_tv    = zeros(1,nWins);
MVI_tv   = zeros(1,nWins);
phi_tv   = zeros(numBands,nWins);

for w = 1:nWins
    idx = (1:winSamps) + (w-1)*stepSamps;
    Xw = EEGDataMat(:,idx);
    Cxw = cov(Xw');
    lam = sort(eig(Cxw),'descend');
    if isempty(lam), lam = eps; end
    FE1_tv(w)  = lam(1)/sum(lam);
    reff_tv(w) = (sum(lam))^2 / sum(lam.^2);

    % band fractions in window
    for b = 1:numBands
        band = Param.Bands(b,:);
        nyq = fs/2;
        W   = band/nyq;
        % Dynamic filter order to avoid filtfilt error
        filtOrder = min(128, floor((size(Xw,2)-1)/3));
        if filtOrder < 1
            Xf = Xw; % skip filtering if too short
        else
            b_filt = fir1(filtOrder, W, 'bandpass');
            Xf = filtfilt(b_filt, 1, double(Xw'))';
        end
        Cb = cov(Xf');
        phi_tv(b,w) = trace(Cb)/trace(Cxw);
    end
    phiMax_w = max(phi_tv(:,w));

    % RID per window
    dR_tv(w) = estimateRID(Xw, fs);

    % MVI per window
    MVI_tv(w) = (FE1_tv(w)/reff_tv(w)) * (1/(1+(dR_tv(w)-1))) * phiMax_w;
end

%% -------------------- Store results --------------------
MSResults.Validity.Global.FE1       = FE1;
MSResults.Validity.Global.reff      = reff;
MSResults.Validity.Global.RID       = dR;
MSResults.Validity.Global.PhiBands  = phi_b;
MSResults.Validity.Global.PhiMax    = phi_max;
MSResults.Validity.Global.MVI       = MVI;
MSResults.Validity.Global.Decision  = ValidityDecision;

MSResults.Validity.TimeVarying.FE1      = FE1_tv;
MSResults.Validity.TimeVarying.reff     = reff_tv;
MSResults.Validity.TimeVarying.RID      = dR_tv;
MSResults.Validity.TimeVarying.PhiBands = phi_tv;
MSResults.Validity.TimeVarying.MVI      = MVI_tv;
MSResults.Validity.TimeVarying.Tvec     = ((1:nWins)-1)*stepMs/1000;

end

%% -------------------- Helper function --------------------
function dR = estimateRID(DataMat, fs)
% Entropy-scaling Rényi Information Dimension (RID) proxy
winMs   = 200; stepMs = 50; nScales = 8;
winSamps  = round(winMs*fs/1000);
stepSamps = round(stepMs*fs/1000);
scales = logspace(-2,0,nScales);

[nch, ns] = size(DataMat);
X = bsxfun(@minus, DataMat, mean(DataMat,2));
X = bsxfun(@rdivide, X, std(X,0,2)+eps);

slopes = [];
for start = 1:stepSamps:(ns-winSamps+1)
    win = X(:, start:start+winSamps-1);
    Hs = zeros(1,nScales);
    for si = 1:nScales
        eps_val = scales(si)*median(std(win,0,2));
        if eps_val <= 0, continue; end
        q = floor(win/eps_val);
        [~,~,ic] = unique(q','rows');
        p = histcounts(ic,1:(max(ic)+1))/length(ic);
        p = p(p>0);
        Hs(si) = -sum(p.*log(p+eps));
    end
    if all(Hs==0), continue; end
    coeffs = polyfit(log(1./scales), Hs, 1);
    slopes(end+1) = coeffs(1);
end
if isempty(slopes)
    dR = 1;
else
    dR = max(1, median(slopes));
end
end
