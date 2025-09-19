function [results] = FnMicrostateValidity_V2(EEGMat, fs)
% FnMicrostateValidity - Compute time-varying EEG features and quantify band relevance for microstates
% Author: Ehsan Eqlimi
% Affiliation: Ghent University
% Email: ehsan.eqlimi@ugent.be
%
% Inputs:
%   EEGMat - EEG data matrix (n_channels x T_samples)
%   fs - Sampling frequency (Hz)
% Outputs:
%   results - Struct containing time-varying features and band relevance:
%       .time - Time points of windows (seconds)
%       .FE1 - Fraction of variance explained by first eigenvalue
%       .r_eff - Effective rank of covariance
%       .RID - Rényi Information Dimension proxy
%       .phi_b - Spectral variance fractions for each band
%       .GFP - Global Field Power time series
%       .GFP_var - GFP variance per window
%       .MVI - Microstate Validity Index
%       .band_names - Names of frequency bands
%       .band_relevance - Relevance scores for each band (higher indicates stronger relation to microstates)
%
% Defaults:
%   window_len_s = 0.2 s (200 ms window)
%   step_size_s = 0.05 s (50 ms step)
%   bands = theta [4-7 Hz], alpha [8-12 Hz], beta [13-30 Hz], gamma [31-45 Hz]
%   num_bins = 50 for RID proxy quantization
%
% Dependencies: Signal Processing Toolbox (butter, filtfilt), Statistics Toolbox (corr)

%% Input Validation
validateattributes(EEGMat, {'numeric'}, {'2d'}, 'FnMicrostateValidity', 'EEGMat');
validateattributes(fs, {'numeric'}, {'scalar', 'positive'}, 'FnMicrostateValidity', 'fs');

%% Default Parameters
window_len_s = 0.2; % Window length in seconds
step_size_s = 0.05; % Step size in seconds
bands = struct('theta', [4, 7], 'alpha', [8, 12], 'beta', [13, 30], 'gamma', [31, 45]);
num_bins = 50; % Number of bins for RID quantization (replaces epsilon_scales)

%% Initialize Parameters
[n, T] = size(EEGMat); % n: channels, T: time samples
window_len = round(window_len_s * fs); % Window length in samples
step_size = round(step_size_s * fs); % Step size in samples
num_windows = floor((T - window_len) / step_size) + 1;
band_names = fieldnames(bands);

%% Preprocess EEG Data
% Center the data (subtract mean across channels)
EEGMat = EEGMat - mean(EEGMat, 1);

%% Initialize Output Structure
results = struct();
results.time = (0:num_windows-1) * step_size / fs; % Time points (seconds)
results.FE1 = zeros(1, num_windows);
results.r_eff = zeros(1, num_windows);
results.RID = zeros(1, num_windows);
results.phi_b = zeros(num_windows, length(band_names));
results.GFP = zeros(1, T);
results.GFP_var = zeros(1, num_windows);
results.MVI = zeros(1, num_windows);
results.band_names = band_names;
results.band_relevance = zeros(1, length(band_names)); % Band relevance

%% Compute Time-Varying Features
band_ranks = zeros(num_windows, length(band_names)); % Store band-specific effective ranks
for w = 1:num_windows
    % Extract window
    idx = (w-1) * step_size + 1 : min((w-1) * step_size + window_len, T);
    X_w = EEGMat(:, idx);
    
    % Normalize data to unit variance per channel
    X_w = X_w ./ (std(X_w, 0, 2) + eps); % Avoid division by zero
    
    % Compute broadband covariance
    C_x = (X_w * X_w') / size(X_w, 2);
    [V, D] = eig(C_x);
    lambda = sort(diag(D), 'descend');
    
    % FE1: Fraction of variance explained by first eigenvalue
    results.FE1(w) = lambda(1) / sum(lambda);
    
    % Effective rank
    results.r_eff(w) = (sum(lambda)^2) / sum(lambda.^2);
    
    % RID proxy: Entropy scaling with adaptive quantization
    entropies = zeros(1, num_bins);
    epsilon_scales = linspace(0.1, 2, num_bins) * std(X_w(:)); % Adaptive scales based on data std
    for e = 1:num_bins
        epsilon = epsilon_scales(e);
        X_quant = floor(X_w / epsilon);
        [~, counts] = unique(X_quant', 'rows');
        probs = histcounts(counts, 'Normalization', 'probability');
        entropies(e) = -sum(probs .* log2(probs + eps));
    end
    % Check for sufficient entropy variation
    if std(entropies) > 1e-3
        p = polyfit(log(1./epsilon_scales), entropies, 1);
        results.RID(w) = max(p(1), 0); % Ensure non-negative RID
    else
        results.RID(w) = NaN; % Flag unreliable estimate
        warning('Window %d: Insufficient entropy variation for RID estimation.', w);
    end
    
    % Spectral variance fractions and band-specific ranks
    phi_b = zeros(1, length(band_names));
    for b = 1:length(band_names)
        band = bands.(band_names{b});
        [filt_b, filt_a] = butter(4, band / (fs/2), 'bandpass');
        X_b = filtfilt(filt_b, filt_a, X_w')';
        C_b = (X_b * X_b') / size(X_b, 2);
        phi_b(b) = trace(C_b) / trace(C_x);
        % Compute band-specific effective rank
        [~, D_b] = eig(C_b);
        lambda_b = sort(diag(D_b), 'descend');
        band_ranks(w, b) = (sum(lambda_b)^2) / sum(lambda_b.^2);
    end
    results.phi_b(w, :) = phi_b;
    
    % GFP variance in window
    GFP_w = sqrt(sum(X_w.^2) / n);
    results.GFP_var(w) = var(GFP_w);
end

% Compute full GFP time series
results.GFP = sqrt(sum(EEGMat.^2) / n);

%% Compute Microstate Validity Index (MVI)
for w = 1:num_windows
    FE1 = results.FE1(w);
    r_eff = results.r_eff(w);
    d_R = results.RID(w);
    if isnan(d_R)
        d_R = 1; % Default to low dimension if RID estimation fails
    end
    max_phi_b = max(results.phi_b(w, :));
    results.MVI(w) = (FE1 / r_eff) * (1 / (1 + (d_R - 1))) * max_phi_b;
end

%% Quantify Band Relevance for Microstates
% Relevance score: mean(phi_b) * corr(phi_b, MVI) * (1/mean(band_rank))
% Higher score indicates stronger relation to microstate dynamics
for b = 1:length(band_names)
    mean_phi_b = mean(results.phi_b(:, b));
    corr_phi_mvi = corr(results.phi_b(:, b), results.MVI', 'Type', 'Pearson', 'Rows', 'complete');
    mean_band_rank = mean(band_ranks(:, b));
    % Ensure positive correlation (negative correlation suggests anti-alignment)
    if isnan(corr_phi_mvi) || corr_phi_mvi < 0
        corr_phi_mvi = 0;
    end
    % Normalize rank contribution (lower rank is better for microstates)
    rank_score = 1 / (1 + mean_band_rank); % Inverse rank, capped to avoid infinity
    results.band_relevance(b) = mean_phi_b * corr_phi_mvi * rank_score;
end
% Normalize relevance scores to sum to 1
results.band_relevance = results.band_relevance / sum(results.band_relevance + eps);

%% Display Band Relevance
fprintf('\n=== Band Relevance for Microstates ===\n');
for b = 1:length(band_names)
    fprintf('%s Band: Relevance Score = %.3f\n', band_names{b}, results.band_relevance(b));
end
[max_score, max_idx] = max(results.band_relevance);
fprintf('Most Relevant Band: %s (Score = %.3f)\n', band_names{max_idx}, max_score);

%% Plot Results
figure('Name', 'Time-Varying EEG Features', 'Position', [100, 100, 1200, 800]);
subplot(3, 2, 1);
plot(results.time, results.FE1, 'LineWidth', 1.5);
title('FE1 Over Time'); xlabel('Time (s)'); ylabel('FE1'); grid on;

subplot(3, 2, 2);
plot(results.time, results.r_eff, 'LineWidth', 1.5);
title('Effective Rank Over Time'); xlabel('Time (s)'); ylabel('r_eff'); grid on;

subplot(3, 2, 3);
plot(results.time, results.RID, 'LineWidth', 1.5);
title('RID Proxy Over Time'); xlabel('Time (s)'); ylabel('RID'); grid on;

subplot(3, 2, 4);
plot(results.time, results.phi_b, 'LineWidth', 1.5);
title('Spectral Variance Fractions'); xlabel('Time (s)'); ylabel('\phi_b');
legend(band_names, 'Location', 'best'); grid on;

subplot(3, 2, 5);
plot(results.time, results.GFP_var, 'LineWidth', 1.5);
title('GFP Variance Over Time'); xlabel('Time (s)'); ylabel('GFP Variance'); grid on;

subplot(3, 2, 6);
plot(results.time, results.MVI, 'LineWidth', 1.5);
title('MVI Over Time'); xlabel('Time (s)'); ylabel('MVI'); grid on;

% Plot GFP time series
figure('Name', 'Global Field Power', 'Position', [100, 100, 600, 400]);
plot((0:T-1)/fs, results.GFP, 'LineWidth', 1.5);
title('GFP Time Series'); xlabel('Time (s)'); ylabel('GFP (\muV)'); grid on;

% Plot Band Relevance
figure('Name', 'Band Relevance for Microstates', 'Position', [100, 100, 600, 400]);
bar(categorical(band_names), results.band_relevance);
title('Band Relevance for Microstates'); xlabel('Frequency Band'); ylabel('Relevance Score');
grid on;

end