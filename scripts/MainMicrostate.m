function MainMicrostate
 %% Convert to FieldTrip structure
            % Specify FT settings
            fieldbox='raw';% Define the fieldbox to be used for analysis
            transform='none';  % Specify no transformation for the data
            EEGData=eeglab2fieldtrip(ALLEEG,fieldbox,transform);  % Convert the EEG data from EEGLAB to FieldTrip format
            EEGData.elec.chanpos=EEGData.elec.elecpos;  % Set the electrode positions
            EEGData.unit='mm';  % Set the unit for the electrode positions to millimeters
            elec=EEGData.elec;  % Extract the electrode structure from the EEGData
            cfg=[];  % Initialize configuration structure
            cfg.elec=elec;  % Assign the manually created 'elec' structure to the configuration
            Layout = ft_prepare_layout(cfg);  % Prepare the layout for plotting the electrode positions
            %% Resample to 512 Hz if the sampling frequency is not already 512 Hz
            if EEGData.fsample~=512
                cfg=[];
                cfg.resamplefs=512;
                EEGData=ft_resampledata(cfg,EEGData);
            end
            %% Compute average reference
            cfg=[];  % Reset the configuration structure
            cfg.reref=Param.DoReRef;  % Indicate that the reference needs to be changed (to average reference)
            cfg.refchannel='all';  % Use all channels for the average reference
            cfg.demean=Param.DoDeMean;  % Demean the data (remove mean from each channel)
            cfg.detrend=Param.DoDeTrend;  % Detrend the data (remove linear trends)
            PreprocessedEEGData=ft_preprocessing(cfg,EEGData);  % Apply the preprocessing steps (average reference, demean, detrend)
            EEGData=[];
            %% Filter data
            cfg=[];% Reset the configuration structure
            cfg.bpfilter=Param.DoBandPassFilter;  % Apply a bandpass filter
            cfg.bpfilttype=Param.BandPassFilterType;  % For example, use a FIR windowed-sinc filter type
            cfg.bpfreq=Param.BandPassFrequencyBand;  % Set the bandpass filter frequency range (1-45 Hz)
            FilteredPreprocessedEEGData=ft_preprocessing(cfg,PreprocessedEEGData);  % Apply the bandpass filter to the preprocessed EEG data
            PreprocessedEEGData=[];
            %% Compute GFP (Global Field Potential)
            cfg=[]; % Initialize configuration structure
            cfg.method=Param.GFPCalcMethod;% Choose the method for computing the global mean field (using amplitude)
            GFP=ft_globalmeanfield(cfg,FilteredPreprocessedEEGData);% Compute the GFP for the EEG data using the specified method
            %% Find peak indices
            MinTFdistance=Param.MinPeakDist*FilteredPreprocessedEEGData.fsample/1000;  % Convert minimum peak distance to samples (based on the sample frequency)
            [~,PeakIDx]=findpeaks(GFP.avg,'MinPeakDistance',MinTFdistance);  % Find the indices of peaks in the GFP data
            % Identify and remove noise peaks by comparing their value to the mean + threshold * standard deviation
            NoisePeakIdx=GFP.avg(PeakIDx)>(mean(GFP.avg)+Param.GFPthresh*std(GFP.avg));
            PeakIdx=PeakIDx(~NoisePeakIdx);  % Remove noisy peaks
            %% Select GFP peaks
            % If you wish to select a predefined number of peaks, you can use the commented-out lines below
            % selection=randperm(length(peakidx));  % Randomly shuffle peak indices
            % gfppeakidx=peakidx(selection(1:nPeaks));  % Select a predefined number of peaks (e.g., nPeaks)
            GFPPeakIdx=PeakIdx;  % Select all remaining peaks
            TempDat=cell2mat(FilteredPreprocessedEEGData.trial);  % Convert the trial data to a matrix format
            GFPData=TempDat(:,GFPPeakIdx);  % Append the selected GFP peaks to the data
            TempDat=[];
            FilteredPreprocessedEEGData=[];
            %% STEP 2: Micro-segmenting EEG Data
            % The goal of this step is to segment the EEG data into microstates,
            % which are defined as short, stable periods of brain activity.
            % We begin by normalizing the data to ensure comparability across channels.

            % Normalise EEG GFP data by dividing each channel's data by the average standard deviation across channels.
            % This step ensures that each channel has the same influence in the clustering algorithm.
            if Param.DoNormalizeGFP
                GFPData=GFPData./mean(std(GFPData,0,2));  % Normalize the GFP data
            end
            %% Use modified k-means algorithm for 1st level clustering
            % We apply the k-means clustering algorithm to the GFP data,
            % trying different values of k (number of clusters) from 4 to 10.
            % The modified k-means algorithm improves performance and reliability over traditional k-means.
            % Run the modified k-means algorithm to find the prototypes (cluster centroids)
            % and assign labels to each data point (which segment of data belongs to which cluster).
            [Prototypes,Labels,Results]=modkmeans(GFPData,Param.KRange,Param.Opts);
            % Note: usable output is is Results.A_all
            % Optionally, save the clustering results if needed for further analysis
            if Param.SaveMSResults
                save([MSResultPath SesseionName  '_ModkmeansOutput.mat' ],'Prototypes','Labels','Results');
            end
            %% Add my custom microstate plot
            % This section can be used for visualizing the results of the microstate analysis.
            % Uncomment the topoplot code to generate topographic plots for each cluster.
            % topoplot
            % Create a copy of the original EEG data to store microstate results
            MSEEG=ALLEEG;
            % Store the results from the modified k-means clustering in the MSEEG structure
            % This includes the clustering results and the settings used for clustering
            MSEEG.microstate.Res=Results;  % Store clustering results in the 'microstate' field
            MSEEG.microstate.algorithm_settings.Nmicrostates=Param.KRange;  % Store the range of microstates tested (K_range)
            % Create a figure to visualize the microstates topographicallyBu
            % The 'MicroPlotTopo' function is used to plot the topographic maps of each microstate
            if Param.PopoUpNonSortedMsFigs
                figure;
                MicroPlotTopo(MSEEG,'plot_range',Param.KRange);  % Visualize topographic maps for the range of microstates
                if Param.DoRdBuColorMap
                    % Check if the brewermap toolbox is available; if not, it will load it
                    ft_hastoolbox('brewermap', 1);  % Ensure that the brewermap toolbox is on the path
                    % Change the colormap for the topographic plots to 'RdBu' and reverse the color order
                    % This colormap is suitable for visualizing EEG microstates with distinct colors.
                    colormap(flipud(brewermap(64, 'RdBu')));  % Set the colormap to RdBu and flip it for better visual contrast
                end
                if Param.SaveAllNonSortedMSFigs
                    saveas(gcf,[MSResultPath SesseionName  '_AllNonSortedMSFigs'],'jpg')
                    close gcf
                end
            else
            end