clc
clear all
%%% This code performs the following:
% 1. Automatically detects the presence of a target 
% 2. Identifies the corresponding range index
% 3. Extracts the unwrapped phase of the slow-time signal
% 4. Performs Matched filtering with the transmitted pilot on the unwrapped phase
% 5. Identifies bit start index
% 6. Extracts bits using FFT over each bit duration of the unwrapped phase
%%%

path = 'C:\Users\Lenovo\OneDrive\Documents\MATLAB_2024a\MovingTargetIMtech\';
binFilePath = 'C:\Users\Lenovo\OneDrive\Documents\MATLAB_2024a\MovingTargetIMtech\Test6\';
if ~isfile([binFilePath 'adc_data.mat'])
    rawDataReaderModified([path '1642.setup.json'],'adc_data','datacube',binFilePath,0);
else
    doNothing = 1;
     % Mat Files already exist.
end
%load([binFilePath 'adc_data.mat']);

% Choose between Phase matched filtering (PMF) or differentiated matched filtering (DMF) for bit start identification
BitStartIdentifyMethod  = "PMF"; % either "PMF" or  "DMF"
BitEstimationMethod     = "FFT"; % either "FFT" or "AR"

load([binFilePath 'adc_data']);
load([binFilePath 'datacube']);

%audioPath           = '..\NSTL_data\WaterExperiments\PoolExpts\audio_filesOFDMpilot0.3NP\';
%[audio,audiofs]     = audioread([audioPath 'SyncPilot1.wav']);
%PilotDuration       = 0.3;
%pilot               = audio(1:round(PilotDuration*audiofs));


mmWaveJsonFileName = [path '1642.mmwave.json'];
mmWaveJSON = jsondecode(fileread(mmWaveJsonFileName));

ConfigParams = mmWaveJSON.mmWaveDevices.rfConfig.rlProfiles.rlProfileCfg_t;

% Radar parameters
c               = 3e8;
fstart          = adcRawData.rfParams.startFreq*1e9;
numADCSamples   = adcRawData.rfParams.numRangeBins;
fs              = adcRawData.rfParams.sampleRate*1e6;
disp(fs);
Slope           = adcRawData.rfParams.freqSlope*1e12; % MHz/us
idleTime        = ConfigParams.idleTimeConst_usec*1e-6;
adc_start_time  = ConfigParams.adcStartTimeConst_usec*1e-6;
f0              = fstart+Slope*adc_start_time;
lambda          = c/f0;
adc_sample_time = numADCSamples/fs;
f1              = fstart+Slope*(adc_start_time+adc_sample_time);
BW              = f1 - f0;
rampEndTime     = ConfigParams.rampEndTime_usec*1e-6;
rangeRes        = c/(2*BW);
Tcri            = rampEndTime + idleTime;
fcrf            = 1/Tcri;
tbar            = (0:numADCSamples-1)/fs;

nch             = radarCube.dim.numRxChan;
numChirps       = radarCube.dim.numChirps;
fst             = (-numChirps/2:numChirps/2-1)/numChirps*fcrf;
f               = (0:numADCSamples-1)/numADCSamples*fs;
% For Positive chirp
if Slope > 0
    ridx            = c*f/(2*Slope);%-c*fs/(2*2*Slope);
else
    ridx            = c*(fs-f)/(2*abs(Slope));
end
nFrames         = radarCube.dim.numFrames;

% Stacking the data across multiple frames in DataAllFrames variable
ch_1F           = [];
DataAllFrames   = [];
tidx            = [];
for kf = 1:nFrames
    tidx            = [tidx; (0:numChirps-1)'*Tcri + (kf-1)*adcRawData.rfParams.framePeriodicity*1e-3];
%     ch_1F           = squeeze(sum(radarCube.data{kf},2));
    ch_1F           = squeeze(radarCube.data{kf}(:,1,:));
    DataAllFrames   = [DataAllFrames; ch_1F];
    
end

backgroundSig = mean(DataAllFrames);
new_data = DataAllFrames - ones(51200,1)*backgroundSig;

nFFT_Range = 256;
nFFT_Doppler = 128;
window_size = 10;
RDM_avg=zeros(numChirps,numADCSamples);
for frameIdx = 1:400                           
    data_points = [];

    %for k=frameIdx:(frameIdx+window_size-1)
        RDM = zeros(nFFT_Doppler,nFFT_Range,nch);

        for rxIdx = 1:nch
            rxData = squeeze(radarCube.data{frameIdx}(:,rxIdx,:));
            %rxData = rxData - mean(rxData,1);
        
            % RDM(:,:,rxIdx) = response(rxData);
            winHann         = hann(numChirps)*ones(1,numADCSamples);
            RDM(:,:,rxIdx) = fftshift(fft(winHann.*rxData,[],1),1);
            %RDM(:,:,rxIdx) = fftshift(fft(rxData,[],1),1);
        end

        RDM_avg = mean(RDM,3);
        RDM_dB = 10*log10(abs(RDM_avg)+eps);
        RDM_linear = 10.^(RDM_dB/10);

        fig = figure('Visible','off');
        imagesc(ridx,((-nFFT_Doppler/2):nFFT_Doppler/2-1)/nFFT_Doppler*fcrf,RDM_linear);hold on;
        reverse_ridx = ridx(end:-1:1);
        ylabel('Doppler Frequency (Hz)');
        xlabel('Range (m)');
        title(['Range-Doppler Map for Frame ' num2str(frameIdx)]);
        colorbar;
        axis xy;

        if ~exist('RDM_images', 'dir')
            mkdir('RDM_images');
        end

        filename = sprintf('RDM_Frame_%03d.png', frameIdx);
        saveas(fig, fullfile('RDM_images', filename));
        close(fig);
        




        % save('rdm_data.mat', 'RDM_dB');
         freqs = ((-nFFT_Doppler/2):nFFT_Doppler/2-1)/nFFT_Doppler*fcrf;


    %[rows,columns] = size(RDM_dB);
    %RDM_flat = reshape(RDM_dB,[],1);

    %[clusterlabels,~] = kmeans(RDM_flat,2);
    %clusteredimage = reshape(clusterlabels,rows,columns);
    %figure;
    %imagesc(ridx, ((-nFFT_Doppler/2):nFFT_Doppler/2-1)/nFFT_Doppler*fcrf, clusteredimage);
    %colorbar;
    %title('K-means Clustering on Range-Doppler Map');
    %xlabel('Range (m)');
    %ylabel('Doppler Frequency (Hz)');
    %axis xy;
    
        cfar2d = phased.CFARDetector2D('GuardBandSize',[3 3],'TrainingBandSize',[5 5],...
        'ProbabilityFalseAlarm',1e-5,'OutputFormat','Detection index');


        cutIdx = [];
        for f_idx = 11:length(freqs)-10
            for r_idx = 11:length(ridx)-10
                cutIdx = [cutIdx,[f_idx;r_idx]];
            end
        end
 
        
        % csvFileName = sprintf('RDM_linear_Frame%d.csv', frameIdx);
        % writematrix(RDM_linear, fullfile(binFilePath, csvFileName));

        
        detections = cfar2d(RDM_linear,cutIdx);

        % imagesc(ridx,((-nFFT_Doppler/2):nFFT_Doppler/2-1)/nFFT_Doppler*fcrf,RDM_linear);hold on;
        % reverse_ridx = ridx(end:-1:1);
        % plot(ridx(detections(2,:)),freqs(detections(1,:)),'r+','MarkerSize',30);
        % ylabel('Doppler Frequency (Hz)');
        % xlabel('Range (m)');
        % title(['Range-Doppler Map with CFAR Detections for Frame ' num2str(frameIdx)]);
        % colorbar;
        % axis xy;
        % pause(1);
        
        %data_points = RDM_dB(:);
        %data_points = [detections(2,:);detections(1,:)]';
    %end

    % doppler_bins = detections(1,:)';
    % range_bins = detections(2,:)';
    % 
    % doppler_vals = freqs(doppler_bins);
    % range_vals = ridx(range_bins);
    % 
    % X = [range_vals, doppler_vals]';
    % 
    % k = 2;
    % GMMModel = fitgmdist(X,k,'RegularizationValue',1e-5);
    % cluster_labels = cluster(GMMModel, X);
    % 
    % figure;
    % gscatter(X(:,1), X(:,2), cluster_labels);  % Range vs Doppler
    % xlabel('Range (m)');
    % ylabel('Doppler (Hz)');
    % title('GMM Clustering on CFAR Detections');
    % grid on;

    

    % numRangeBins = 20;
    % numDopplerBins = 20; 
    % 
    % % Compute grid boundaries
    % rangeEdges = linspace(min(ridx), max(ridx), numRangeBins + 1);
    % dopplerEdges = linspace(min(freqs), max(freqs), numDopplerBins + 1);
    % 
    % % Initialize grid for summing RDM_linear intensities
    % gridIntensity = zeros(numDopplerBins, numRangeBins);
    % 
    % % Map RDM_linear values to grid cells
    % [rows, cols] = size(RDM_linear);
    % for r = 1:rows
    %     for c = 1:cols
    %         rangeVal = ridx(c); % Range value for column c
    %         dopplerVal = freqs(r); % Doppler value for row r
    % 
    %         % Find corresponding grid cell
    %         rangeBin = find(rangeEdges(1:end-1) <= rangeVal & rangeVal < rangeEdges(2:end), 1);
    %         dopplerBin = find(dopplerEdges(1:end-1) <= dopplerVal & dopplerVal < dopplerEdges(2:end), 1);
    % 
    %         % Add RDM_linear intensity to grid cell
    %         if ~isempty(rangeBin) && ~isempty(dopplerBin)
    %             gridIntensity(dopplerBin, rangeBin) = gridIntensity(dopplerBin, rangeBin) + RDM_linear(r, c);
    %         end
    %     end
    % end
    % 
    % % Normalize grid intensity by the number of RDM_linear pixels per cell
    % rangeBinWidth = diff(rangeEdges(1:2));
    % dopplerBinWidth = diff(dopplerEdges(1:2));
    % pixelsPerCell = (rangeBinWidth / (max(ridx) - min(ridx)) * cols) * ...
    %             (dopplerBinWidth / (max(freqs) - min(freqs)) * rows);
    % gridIntensity = gridIntensity / pixelsPerCell; % Average intensity per cell
    % 
    % % Threshold for identifying significant grid cells
    % intensityThreshold = prctile(gridIntensity(:), 95); % Top 5% of intensities
    % significantCells = gridIntensity >= intensityThreshold;
    % 
    % % Label connected components (clusters) in the grid
    % [labeledGrid, numClusters] = bwlabel(significantCells, 4); % 4-connectivity
    % 
    % % Extract cluster centroids for visualization
    % clusterCentroids = zeros(numClusters, 2); % [range, doppler]
    % for c = 1:numClusters
    %     [dopplerBins, rangeBins] = find(labeledGrid == c);
    %     if ~isempty(dopplerBins)
    %         % Compute centroid in range-doppler space
    %         rangeCentroid = mean(rangeEdges(rangeBins));
    %         dopplerCentroid = mean(dopplerEdges(dopplerBins));
    %         clusterCentroids(c, :) = [rangeCentroid, dopplerCentroid];
    %     end
    % end
    % 
    % % Visualize the Range-Doppler Map with Grid-Based Clustering
    % figure;
    % imagesc(ridx, freqs, RDM_dB); hold on;
    % colormap('jet');
    % colorbar;
    % xlabel('Range (m)');
    % ylabel('Doppler Frequency (Hz)');
    % title(['Grid-Based Clustering on RDM for Frame ' num2str(frameIdx)]);
    % axis xy;
    % 
    % % Plot cluster centroids
    % colors = {'rx', 'go', 'b+', 'm*', 'cs'}; % Colors for different clusters
    % for c = 1:numClusters
    %     if ~isnan(clusterCentroids(c, 1))
    %         scatter(clusterCentroids(c, 1), clusterCentroids(c, 2), 50, ...
    %             colors{mod(c-1, length(colors))+1}, 'LineWidth', 2, ...
    %             'DisplayName', ['Cluster ' num2str(c)]);
    %     end
    % end
    % legend('show');
    % hold off;
    % pause(1); % Pause to view the plot

    

    % % IQR Method to calculate the mask
    % Q1 = quantile(data_points(:), 0.25);
    % Q3 = quantile(data_points(:), 0.75);
    % IQR_value = Q3 - Q1;
    % lower_bound = Q1 - 1.5 * IQR_value;
    % upper_bound = Q3 + 1.5 * IQR_value;
    % 
    % % Outlier removal by removing points outside the mask
    % mask = (data_points >= lower_bound) & (data_points <= upper_bound);
    % data_points(~mask) = NaN;  % Mark outliers as NaN
    % filtered_data = data_points(~isnan(data_points));
    % 
    % % Standardization
    % mu = mean(filtered_data);
    % sigma = std(filtered_data);
    % normalized_data_points = (filtered_data - mu) / sigma;

    [rows, cols] = size(RDM_linear);
    [rowIdx, colIdx] = ind2sub([rows, cols], 1:numel(RDM_linear));

    % KMeans Clustering
    % k = 2;
    % [idx,C] = kmeans(normalized_data_points,k);
    % figure;
    % imagesc(ridx, freqs, RDM_linear); hold on;
    % xlabel('Range (m)');
    % ylabel('Doppler Frequency (Hz)');
    % title(['K-means Clustering on RDM for frame ' num2str(frameIdx)]);
    % axis xy;
    % scatter(data_points(idx == 1,1), data_points(idx == 1,2), 100, 'rx', 'LineWidth', 2);
    % scatter(data_points(idx == 2,1), data_points(idx == 2,2), 100, 'go', 'LineWidth', 2);
    % legend('Target 1', 'Target 2');
    % pause(1)
    % hold off;

    % DBSCAN Clustering
    % epsilon = 0.5;
    % minPts = 3;
    % idx = dbscan(normalized_data_points, epsilon, minPts);
    % figure;
    % imagesc(ridx, freqs, RDM_dB); hold on;
    % xlabel('Range (m)');
    % ylabel('Doppler Frequency (Hz)');
    % title(['DBSCAN Clustering on CFAR Detections for frame ' num2str(frameIdx)]);
    % axis xy;
    % scatter(data_points(idx == 1,1), data_points(idx == 1,2), 100, 'rx', 'LineWidth', 2);
    % scatter(data_points(idx == 2,1), data_points(idx == 2,2), 100, 'go', 'LineWidth', 2);
    % legend('Cluster 1', 'Cluster 2', 'Noise');
    % pause(1);
    % hold off;

    % GMM Clustering
    % gmmModel = fitgmdist(normalized_data_points,2);
    % idx = cluster(gmmModel,normalized_data_points);
    % figure;
    % imagesc(ridx,freqs,RDM_dB);hold on;
    % xlabel('Range (m)');
    % ylabel('Doppler Frequency (Hz)');
    % title(['GMM CLustering on CFAR Detections for frame ' num2str(frameIdx)]);
    % colorbar;
    % axis xy;
    % scatter(data_points(idx == 1), data_points(idx == 1), 100, 'rx', 'LineWidth', 2);
    % scatter(data_points(idx == 2), data_points(idx == 2), 100, 'go', 'LineWidth', 2);
    % legend('Target 1', 'Target 2', 'GMM Centroids');
    % pause(1)
    % hold off;

    % Spectral Clustering
    % affinity_matrix = squareform(pdist(normalized_data_points));
    % [idx,~] = spectralcluster(affinity_matrix,2);
    % figure;
    % imagesc(ridx,freqs,RDM_dB);hold on;
    % xlabel('Range (m)');
    % ylabel('Doppler Frequency (Hz)');
    % title(['Spectral Clustering on CFAR Detections for frame ' num2str(frameIdx)]);
    % colorbar;
    % axis xy;
    % scatter(data_points(idx==1,1),data_points(idx==1,2),100,'rx','LineWidth',2);
    % scatter(data_points(idx==2,1),data_points(idx==2,2),100,'go','LineWidth',2);
    % legend('Cluster 1','Cluster 2');
    % pause(1);
    % hold off;

    % Hierarchical KMeans Clustering
    %k = 2;
    %[idx, ~] = kmeans(normalized_data_points, k);
    %SVMModel = fitcsvm(normalized_data_points, idx, 'KernelFunction', 'linear');
    %predicted_labels = predict(SVMModel, normalized_data_points);
    % Plot Range-Doppler Map
    %figure;
    %imagesc(ridx, freqs, RDM_dB); hold on;
    %xlabel('Range (m)');
    %ylabel('Doppler Frequency (Hz)');
    %title('SVM Classification on CFAR Detections');
    %axis xy;
    %scatter(data_points(predicted_labels == 1,1), data_points(predicted_labels == 1,2), 100, 'rx', 'LineWidth', 2);
    %scatter(data_points(predicted_labels == 2,1), data_points(predicted_labels == 2,2), 100, 'go', 'LineWidth', 2);
    %legend('Class 1', 'Class 2');
    %pause(1);
    %hold off;
end