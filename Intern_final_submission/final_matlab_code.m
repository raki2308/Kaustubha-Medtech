% Initialize a file for storing heartbeat information
heartbeatFile = fopen('heartbeat_data.txt', 'w');

previousData = [];  % Initialize variable to store previous data

while true
    % Load ECG data from the CSV file (replace 'Real_FECG_1.csv' with the actual file name)
    currentData = csvread('Real_FECG_2.csv');

    % Check if the data has changed
    if ~isequal(currentData, previousData)
        % Update previous data
        previousData = currentData;
 
    % Load ECG data from the CSV file (replace 'Real_FECG_1.csv' with the actual file name)
    data = csvread('Real_FECG_2.csv');
    ann_data  = csvread('Real_FECG_2_ann.csv');
    
    % Extract the two columns from the data
    column1 = data(:, 1);
    column2 = data(:, 2);
    
    % Adaptive Baseline Correction using Moving Average
    window_size = 100; % Adjust the window size as needed
    
    baseline_corrected_signal1 = zeros(size(column1));
    baseline_corrected_signal2 = zeros(size(column2));
    
    for i = 1:length(column1)
        start_index = max(1, i - window_size);
        end_index = min(length(column1), i + window_size);
        
        % Calculate the local mean (baseline) within the window for each column
        local_mean1 = mean(column1(start_index:end_index));
        local_mean2 = mean(column2(start_index:end_index));
        
        % Subtract the local mean from the signal to correct the baseline
        baseline_corrected_signal1(i) = column1(i) - local_mean1;
        baseline_corrected_signal2(i) = column2(i) - local_mean2;
    end
    
    % Sample rate (you may need to adjust this based on your signal)
    fs = 1500;
    
    % Define parameters
    order = 4; % Filter order
    center_frequency = 50 / (fs/2); % Center frequency of the notch (normalized to Nyquist frequency)
    bandwidth = 2 / (fs/2); % Bandwidth of the notch (normalized to Nyquist frequency)
    
    % Create a bandstop (notch) Butterworth filter
    [bs1_b, bs1_a] = butter(order, [center_frequency - bandwidth/2, center_frequency + bandwidth/2], 'stop');
    
    % Apply the filter to each column of the signal
    filtered_column1bs1 = filter(bs1_b, bs1_a, baseline_corrected_signal1);
    filtered_column2bs1 = filter(bs1_b, bs1_a, baseline_corrected_signal2);
    
    % Define parameters
    order = 4; % Filter order
    center_frequency = 60 / (fs/2); % Center frequency of the notch (normalized to Nyquist frequency)
    bandwidth = 2 / (fs/2); % Bandwidth of the notch (normalized to Nyquist frequency)
    
    % Create a bandstop (notch) Butterworth filter
    [bs2_b, bs2_a] = butter(order, [center_frequency - bandwidth/2, center_frequency + bandwidth/2], 'stop');
    
    % Apply the filter to each column of the signal
    filtered_column1bs2 = filter(bs2_b, bs2_a, filtered_column1bs1);
    filtered_column2bs2 = filter(bs2_b, bs2_a, filtered_column2bs1);
    
    % Define parameters
    order = 4; % Filter order
    cutoff_frequency = 0.33 / (fs/2); % Cutoff frequency (normalized to Nyquist frequency)
    
    % Create a high-pass Butterworth filter
    [h_b, h_a] = butter(order, cutoff_frequency, 'high');
    
    % Apply the filter to each column of the signal
    filtered_column1h = filter(h_b, h_a, filtered_column1bs2);
    filtered_column2h = filter(h_b, h_a, filtered_column2bs2);
    
    % Define parameters
    order = 4; % Filter order
    cutoff_frequency = 150 / (fs/2); % Cutoff frequency (normalized to Nyquist frequency)
    
    % Create a low-pass Butterworth filter
    [l_b, l_a] = butter(order, cutoff_frequency, 'low');
    
    % Apply the filter to each column of the signal
    filtered_column1l = filter(l_b, l_a, filtered_column1h);
    filtered_column2l = filter(l_b, l_a, filtered_column2h);

    % Extract the ECG signals from both columns
    ecg_signal1 = filtered_column1l;
    ecg_signal2 = filtered_column2l;

    csvwrite("Book2.csv",ecg_signal1);
    
        % Apply a moving average filter to smooth the ECG signal
        smoothed_ecg1 = movmean(ecg_signal1, 15);
        
        % Plot the original and smoothed signals
        figure;
        plot(ecg_signal1);
        hold on;
        plot(smoothed_ecg1);
        title('Original and Smoothed ECG Signals');
        xlabel('Sample Index');
        ylabel('Amplitude');
        legend('Original Signal', 'Smoothed Signal');
        
        % Apply a moving average filter to smooth the ECG signal
        smoothed_ecg2 = movmean(ecg_signal2, 15);
        
        % Plot the original and smoothed signals
        figure;
        plot(ecg_signal2);
        hold on;
        plot(smoothed_ecg2);
        title('Original and Smoothed ECG Signals');
        xlabel('Sample Index');
        ylabel('Amplitude');
        legend('Original Signal', 'Smoothed Signal');
    
    % Calculate a dynamic threshold as a percentage of the maximum amplitude
    threshold_percentage = 70;  % Adjust this percentage as needed
    threshold = threshold_percentage / 100 * max(abs(ecg_signal1));
        
    % Find R-peaks in the ECG signal
    [peaks1, peak_indices1] = findpeaks(-ecg_signal1, 'MinPeakHeight', threshold, 'MinPeakDistance', fs/100);
    [peaks2, peak_indices2] = findpeaks(-ecg_signal2, 'MinPeakHeight', threshold, 'MinPeakDistance', fs/100);

    csvwrite('maternal.csv', peak_indices1)
    csvwrite('maternal2.csv', peak_indices2)

    % Calculate heart rate for both columns
    heart_rate1 = calculateHeartRateFromRPeaks(peak_indices1, fs,ecg_signal1);
    heart_rate2 = calculateHeartRateFromRPeaks(peak_indices2, fs,ecg_signal2);
    
    % Display heart rates
    disp(['Heart Rate (Column 1): ' num2str(heart_rate1) ' beats per minute']);
    disp(['Heart Rate (Column 2): ' num2str(heart_rate2) ' beats per minute']);

% Load mixed ECG data from CSV file
mixedECGData = csvread('Real_FECG_2.csv');
ann_data = csvread('Real_FECG_2_ann.csv');

% Define parameters
samplingRate = 2500;      % Replace with the actual sampling rate of your signal
notchFrequency = 50;      % Maternal ECG frequency (Hz)
alpha = 0.1;              % Adaptation parameter

% Apply low-pass filter (cutoff at 150 Hz)
cutoffLowPass = 150;
[bLowPass, aLowPass] = butter(4, cutoffLowPass / (samplingRate / 2), 'low');
mixedECGData = filtfilt(bLowPass, aLowPass, mixedECGData);

% Apply high-pass filter (cutoff at 0.33 Hz)
cutoffHighPass = 0.33;
[bHighPass, aHighPass] = butter(4, cutoffHighPass / (samplingRate / 2), 'high');
mixedECGData = filtfilt(bHighPass, aHighPass, mixedECGData);

    window_size = 100;
    
    baseline_corrected_signal = zeros(size(mixedECGData));
    
    for i = 1:length(mixedECGData)
        start_index = max(1, i - window_size);
        end_index = min(length(mixedECGData), i + window_size);
        
        % Calculate the local mean (baseline) within the window for each column
        local_mean = mean(mixedECGData(start_index:end_index));
        
        % Subtract the local mean from the signal to correct the baseline
        baseline_corrected_signal(i) = mixedECGData(i) - local_mean;
    end
    
% Baseline correction
mixedECGData = baseline_corrected_signal;

% Initialize variables
maternalEstimate = 0;

% Process the mixed ECG signal
fetalECGData = zeros(size(mixedECGData));

% Notch filter for maternal ECG estimation
for i = 1:length(mixedECGData)
    maternalEstimate = alpha * mixedECGData(i) + (1 - alpha) * maternalEstimate;

    % Subtract maternal estimate from the mixed signal to get fetal ECG
    fetalECGData(i) = mixedECGData(i) - maternalEstimate;
end

csvwrite('fetal.csv', fetalECGData);

    % Load ECG data from the CSV file (replace 'Real_FECG_1.csv' with the actual file name)
    data = csvread('fetal.csv');
    maternal_data = csvread('maternal.csv');
    ann_data  = csvread('Real_FECG_2_ann.csv');
    
    % Extract the two columns from the data
    column1 = data(:, 1);
    column2 = data(:, 2);
    
    % Adaptive Baseline Correction using Moving Average
    window_size = 100; % Adjust the window size as needed
    
    baseline_corrected_signal1 = zeros(size(column1));
    baseline_corrected_signal2 = zeros(size(column2));
    
    for i = 1:length(column1)
        start_index = max(1, i - window_size);
        end_index = min(length(column1), i + window_size);
        
        % Calculate the local mean (baseline) within the window for each column
        local_mean1 = mean(column1(start_index:end_index));
        local_mean2 = mean(column2(start_index:end_index));
        
        % Subtract the local mean from the signal to correct the baseline
        baseline_corrected_signal1(i) = column1(i) - local_mean1;
        baseline_corrected_signal2(i) = column2(i) - local_mean2;
    end
    
    % Sample rate (you may need to adjust this based on your signal)
    fs = 1500;
    k= 20;

    % Define parameters
    order = 4; % Filter order
    center_frequency = 50 / (fs/2); % Center frequency of the notch (normalized to Nyquist frequency)
    bandwidth = 2 / (fs/2); % Bandwidth of the notch (normalized to Nyquist frequency)
    
    % Create a bandstop (notch) Butterworth filter
    [bs1_b, bs1_a] = butter(order, [center_frequency - bandwidth/2, center_frequency + bandwidth/2], 'stop');
    
    % Apply the filter to each column of the signal
    filtered_column1bs1 = filter(bs1_b, bs1_a, baseline_corrected_signal1);
    filtered_column2bs1 = filter(bs1_b, bs1_a, baseline_corrected_signal2);

    % Define parameters
    order = 4; % Filter order
    center_frequency = 60 / (fs/2); % Center frequency of the notch (normalized to Nyquist frequency)
    bandwidth = 2 / (fs/2); % Bandwidth of the notch (normalized to Nyquist frequency)
    
    % Create a bandstop (notch) Butterworth filter
    [bs2_b, bs2_a] = butter(order, [center_frequency - bandwidth/2, center_frequency + bandwidth/2], 'stop');
    
    % Apply the filter to each column of the signal
    filtered_column1bs2 = filter(bs2_b, bs2_a, filtered_column1bs1);
    filtered_column2bs2 = filter(bs2_b, bs2_a, filtered_column2bs1);

    % Define parameters
    order = 4; % Filter order
    cutoff_frequency = 0.33 / (fs/2); % Cutoff frequency (normalized to Nyquist frequency)
    
    % Create a high-pass Butterworth filter
    [h_b, h_a] = butter(order, cutoff_frequency, 'high');
    
    % Apply the filter to each column of the signal
    filtered_column1h = filter(h_b, h_a, filtered_column1bs2);
    filtered_column2h = filter(h_b, h_a, filtered_column2bs2);
    
    % Define parameters
    order = 4; % Filter order
    cutoff_frequency = 150 / (fs/2); % Cutoff frequency (normalized to Nyquist frequency)
    
    % Create a low-pass Butterworth filter
    [l_b, l_a] = butter(order, cutoff_frequency, 'low');
    
    % Apply the filter to each column of the signal
    filtered_column1l = filter(l_b, l_a, filtered_column1h);
    filtered_column2l = filter(l_b, l_a, filtered_column2h);

    % Extract the ECG signals from both columns
    ecg_signal1 = filtered_column1l;
    ecg_signal2 = filtered_column2l;
    
        % Apply a moving average filter to smooth the ECG signal
        smoothed_ecg1 = movmean(ecg_signal1, 15);
        
        % Apply a moving average filter to smooth the ECG signal
        smoothed_ecg2 = movmean(ecg_signal2, 15);

    % Extract points from maternal data
    exclude_points = maternal_data;
    
    % Range around exclude_points where peaks should be excluded
    exclude_range = fs/k;
    
    % Create a mask to exclude peaks around specified points
    mask = true(size(data));
    for point = exclude_points'
        % Ensure the point is within the valid range
        valid_range = (point - exclude_range <= (1:length(data))) & ((1:length(data)) <= point + exclude_range);
        mask(valid_range) = false;
    end
    
    % Apply the mask to the signal
    masked_signal = data .* mask;

    % Ensure that masked_signal is a column vector
    masked_signal = masked_signal(:);

    % Calculate a dynamic threshold as a percentage of the maximum amplitude
    threshold_percentage1 = 20;  % Adjust this percentage as needed
    threshold1 = threshold_percentage1 / 100 * max(abs(ecg_signal1));
        
    % Find R-peaks in the ECG signal
    [peaks1a, peak_indices1a] = findpeaks(-masked_signal, 'MinPeakHeight', threshold1, 'MinPeakDistance', fs/k);
    [peaks2a, peak_indices2a] = findpeaks(-ecg_signal2, 'MinPeakHeight', threshold1, 'MinPeakDistance', fs/k);

    % Calculate a dynamic threshold as a percentage of the maximum amplitude
    %threshold_percentage2 = 30;  % Adjust this percentage as needed
    %threshold3 = threshold_percentage3 / 100 * max(abs(ecg_signal1))
        
    % Find R-peaks in the ECG signal
    [peaks1i, peak_indices1i] = findpeaks(masked_signal, 'MinPeakHeight', threshold1, 'MinPeakDistance', fs/k);
    [peaks2i, peak_indices2i] = findpeaks(ecg_signal2, 'MinPeakHeight', threshold1, 'MinPeakDistance', fs/k);
    
    % Determine R-peaks based on distance from baseline
    if mean(peaks1i) - mean(peaks1a) > mean(peaks1a) - mean(peaks1i)
        peaks1 = peak_indices1i;
    else
        peaks1 = peak_indices1a;
    end
    
    if mean(peaks2i) - mean(peaks2a) > mean(peaks2a) - mean(peaks2i)
        peaks2 = peak_indices2i;
    else
        peaks2 = peak_indices2a;
    end
    
    thresholds = [threshold1, -threshold1];
    
        % Plot the ECG signal with detected R-peaks
        figure;
        plot(ecg_signal1);
        hold on;
        plot(peaks1, ecg_signal1(peaks1), 'r*', 'MarkerSize', 10);
        hold on;
        for i = 1:length(ann_data)
            x_value = ann_data(i);
            plot([x_value, x_value], ylim, 'k--');  % Draw vertical lines in black dashed style
        end
        hold on;
        for i = 1:length(thresholds)
            y_value = thresholds(i);
            plot(xlim, [y_value, y_value], 'k--');  % Draw horizontal lines in black dashed style
        end
        title('ECG Signal with Detected R-peaks');
        xlabel('Time (s)');
        ylabel('Amplitude');

    % Calculate distances between consecutive points
    distances = diff(peaks1);

    % Calculate the average distance
    avgDistance = mean(distances);

    % Find pairs of points with distances greater than the average
    idxToInsert = find(distances > 1.5 * avgDistance);

    newpoints=zeros(1,length(peaks1)+length(idxToInsert));
    k=0;
    i=1;
    while i+k<=length(newpoints) && i<=length(peaks1)
        if k+1<=length(idxToInsert) && i==idxToInsert(k+1)+1
            newpoints(i+k)=floor(mean(peaks1(i-1:i)));
            k=k+1;
        else
            newpoints(i+k)=peaks1(i);
            i=i+1;
        end
    end
    
    peaks1 = newpoints;

    % Calculate distances between consecutive points
    distances2 = diff(peaks2);

    % Calculate the average distance
    avgDistance2 = mean(distances2);

    % Find pairs of points with distances greater than the average
    idxToInsert2 = find(distances2 > avgDistance2);

    newpoints2=zeros(1,length(peaks2)+length(idxToInsert2));
    k=0;
    i=1;
    while i+k<=length(newpoints2) && i<=length(peaks2)
        if k+1<=length(idxToInsert2) && i==idxToInsert2(k+1)+1
            newpoints2(i+k)=floor(mean(peaks2(i-1:i)));
            k=k+1;
        else
            newpoints2(i+k)=peaks2(i);
            i=i+1;
        end
    end
    
    peaks2 = newpoints2;

    % Function to calculate heart rate
    %calculateHeartRate = @(ecg_signal) calculateHeartRateFromECG(ecg_signal, fs, threshold);

    % Calculate heart rate for both columns
    heart_rate1 = calculateHeartRatefromRPeaks(peaks1, fs,ecg_signal1);
    heart_rate2 = calculateHeartRatefromRPeaks(peaks2, fs,ecg_signal2);
    
    % Display heart rates
    disp(['Heart Rate (Column 1): ' num2str(heart_rate1) ' beats per minute']);
    disp(['Heart Rate (Column 2): ' num2str(heart_rate2) ' beats per minute']);
    
        % Plot the ECG signal with detected R-peaks
        figure;
        plot(ecg_signal1);
        hold on;
        plot(peaks1, ecg_signal1(peaks1), 'r*', 'MarkerSize', 10);
        hold on;
        for i = 1:length(ann_data)
            x_value = ann_data(i);
            plot([x_value, x_value], ylim, 'k--');  % Draw vertical lines in black dashed style
        end
        hold on;
        for i = 1:length(thresholds)
            y_value = thresholds(i);
            plot(xlim, [y_value, y_value], 'k--');  % Draw horizontal lines in black dashed style
        end
        title('ECG Signal with Detected R-peaks');
        xlabel('Time (s)');
        ylabel('Amplitude');
    end
    
        % Pause for 10 seconds before the next iteration
        pause(10);
end

% Close the heartbeat file when done
fclose(heartbeatFile);

% Function to calculate heart rate from R-peaks
    function heart_rate = calculateHeartRatefromRPeaks(peaks, fs, ecg_signal)
        % Calculate the time between consecutive R-peaks to get the heart rate
        RR_intervals = diff(peaks) / fs;
        heart_rate = 60 / mean(RR_intervals);
    end

% Function to calculate heart rate from R-peaks
    function heart_rate = calculateHeartRateFromRPeaks(peak_indices, fs, ecg_signal)
        % Calculate the time between consecutive R-peaks to get the heart rate
        RR_intervals = diff(peak_indices) / fs;
        heart_rate = 60 / mean(RR_intervals);
        
        % Plot the ECG signal with detected R-peaks
        figure;
        plot((0:length(ecg_signal)-1)/fs, ecg_signal);
        hold on;
        plot(peak_indices/fs, ecg_signal(peak_indices), 'r*', 'MarkerSize', 10);
        title('ECG Signal with Detected R-peaks');
        xlabel('Time (s)');
        ylabel('Amplitude');
        
        % Calculate the duration of the ECG signal
        signal_duration = length(ecg_signal) / fs;

        %time = 20;
        % Set x-axis limits
        xlim([0, signal_duration]);
    
        legend('ECG Signal', 'R-peaks');
      end
