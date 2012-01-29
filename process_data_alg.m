function process_data_alg( file_without_extension )
%PROCESS_DATA_ALG Converts raw data into input data for FIT_GAUSS_MIXTURES
%   We find the .MAT file containing the TimeStamps and DataPoints, and
%   extract the desired dimensions from it. Currently, we use the
%   TimeStamps, Peaks for each wire, and PC's 1-3 of each wire. Widths need
%   to be implemented in way that doesn't break function.

fwe = file_without_extension;
pathToFoo = 'C:\Users\Alexander\Desktop\Aaron Fryman\Data'; %'F:\AaronF\Data\Algorithm\'
foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
load(foo);

% figure; plot(dp(:,1:100))

N = 3;
newDP = dp';
nrows = size(newDP, 1);
score = zeros(nrows, 128);
peak_values = zeros(nrows, 4);
valleys = zeros(nrows, 4);
prin_comps = zeros(12, nrows);
bad_wires = [0 0 0 0];

widths = zeros(nrows,4);
nonlinear_energies = zeros(nrows,4);
peak_start = zeros(nrows,4);
peak_end = zeros(nrows,4);
delay_width = zeros(nrows,4);
delay_start = zeros(nrows,4);
delay_end = zeros(nrows,4);
valley_width = zeros(nrows,4);
valley_start = zeros(nrows,4);
valley_end = zeros(nrows,4);

warning off dg_spikewidth:novalley
warning off dg_spikewidth:novalley2
warning off dg_spikewidth:earlypeak

for wire = 1:4
    disp('WIRE STARTED');
    startIdx = (wire - 1)*32 + 1;
    endIdx = 32*wire;
    wire_data = newDP(:, startIdx:endIdx);
    if mean(wire_data) ~= zeros(1, 32)
        [~, score(:, startIdx:endIdx)] = princomp(wire_data);
        peak_values(:, wire) = max(wire_data, [], 2);
        valleys(:,wire) = min(wire_data, [], 2);
        prin_comps((wire - 1)*N + 1:wire*N, :) = score(:, startIdx:startIdx + N - 1)';
        
        for i=1:nrows
            if mod(i,10000)==0
                if strcmp(progress,'..........')==1
                    progress = '.';
                else
                    progress = strcat(progress,'.');
                end
                disp(progress);
            end
            [peak_width_vals, delay_vals, valley_vals] = dg_spikewidth(wire_data(i,:), 'allpoints');
            widths(i,wire) = peak_width_vals(1);
            t_0 = round(peak_width_vals(2));
            t_1 = round(peak_width_vals(3));
            
            peak_start(i,wire) = t_0;
            peak_end(i,wire) = t_1;
            delay_width(i,wire) = delay_vals(1);
            delay_start(i,wire) = delay_vals(2);
            delay_end(i,wire) = delay_vals(3);
            valley_width(i,wire) = valley_vals(1);
            valley_start(i,wire) = valley_vals(2);
            valley_end(i,wire) = valley_vals(3);
            
            if isnan(widths(i,wire)) == 0
                squared_voltage = wire_data(i,t_0:t_1).*wire_data(i,t_0:t_1);
                if t_0 ~= 1 && t_1 ~= 32
                    nonlinear_energies(i,wire) = sum(squared_voltage...
                        - wire_data(i,t_0-1:t_1-1).*wire_data(i,t_0+1:t_1+1))/widths(i,wire);
                elseif t_1 ~= 32
                    nonlinear_energies(i,wire) = sum(squared_voltage...
                        - [wire_data(i,t_0), wire_data(i,t_0:t_1-1)].*wire_data(i,t_0+1:t_1 + 1))/widths(i,wire);
                elseif t_0 ~= 1
                    nonlinear_energies(i,wire) = sum(squared_voltage...
                        - wire_data(i,t_0-1:t_1-1).*[wire_data(i,t_0+1:t_1), wire_data(i,t_1)])/widths(i,wire);
                else
                    nonlinear_energies(i,wire) = sum(squared_voltage...
                        - [wire_data(i,t_0), wire_data(i,t_0:t_1-1)].*[wire_data(i,t_0+1:t_1), wire_data(i,t_1)])/widths(i,wire);
                end
            else
                nonlinear_energies(i,wire) = 0;
            end
        end
    else
        bad_wires(wire) = 1;
    end
end
    
    % smoothwave_max_wire = max(smoothwave);
    % conv_noise_rating = 1./(1+exp(5*10^-4*(-1*smoothwave_max_wire+3*10^4)));
    % weight_coeffs = [1];
    % weight_terms = [conv_noise_rating];
    % weights = weight_coeffs.*weight_terms;
    
    if any(bad_wires)
        for badIdx = find(bad_wires)
            startIdx = (badIdx - 1)*3 + 1;
            endIdx = 3*badIdx;
            peak_values(:, badIdx) = [];
%             valleys(:, badIdx) = [];
            prin_comps(startIdx:endIdx, :) = [];
        end
    end
    
    input = [TimeStamps; peak_values'; prin_comps];
    
    conv_noise_threshold = 0.1;
    energy_noise_threshold = 500;
    
    baseline_noise_idx = logical((conv_filter_rating < conv_noise_threshold) + (max(nonlinear_energies, [], 2)' > energy_noise_threshold));
    noise_feature_data = input(:,baseline_noise_idx);
    cluster_feature_data = input(:,baseline_noise_idx==0);
    
    input_weights = conv_filter_rating(baseline_noise_idx==0);
    reconstructed_noise_rating = [input_weights conv_filter_rating(baseline_noise_idx)];
    reconstructed_input = [cluster_feature_data noise_feature_data];
    reconstructed_waveform_data = [waveform_data(baseline_noise_idx==0, :); waveform_data(baseline_noise_idx, :)];
    reconstructed_noise_feature_data = noise_feature_data;
    reconstructed_cluster_feature_data = cluster_feature_data;
    
    pathToFoo = '~/UROP/Data/Algorithm_Inputs'; %'C:\Users\Alexander\Desktop\Aaron Fryman\Data\Algorithm_Inputs'; %'F:\AaronF\Data\Algorithm\'
    foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
    save(foo, 'reconstructed_noise_rating', 'reconstructed_input', 'reconstructed_waveform_data', 'reconstructed_noise_feature_data', 'reconstructed_cluster_feature_data');
    
end

%         t_0 = zeros(nrows,1); % Were used before for computing widths
%         t_1 = zeros(nrows,1);
%         for i=1:nrows
%             waveform = wire_data(i,:) - mean(wire_data(i,:));
%             abswave = abs(waveform');
%             smoothwave(wire,i) = max(conv(abswave, hw, 'valid'));
%
%             [local_min_idx,~,~] = dg_findFlattops(wire_data(i,:)*(-1));
%             peak_idx = peak_indices(i);
%             relative_idx = sort([1; local_min_idx; peak_idx; 32]);
%             if peak_idx < 2 || peak_idx > 31
%                 t_0(i) = 0;
%                 t_1(i) = 0;
%                 velocities(i,1,wire) = 0;
%                 velocities(i,2,wire) = 0;
%                 energies(i,wire) = 0;
%             nonlinear_energies(i,wire) = 0;

%                 t_0 = relative_idx(find(relative_idx == peak_idx)-1);
%                 t_1 = relative_idx(find(relative_idx == peak_idx)+1);
%                 widths(i,wire) = t_1-t_0;
%                 velocities(i,1,wire) = (peak_values(i) - wire_data(i,t_0))/(peak_idx - t_0);
%                 velocities(i,2,wire) = (wire_data(i,t_1) - peak_values(i))/(t_1 - peak_idx);
%                 energies(i,wire) = sum(wire_data(i,t_0:t_1).*wire_data(i,t_0:t_1))/widths(i,wire);

%             [~,t_0(i)] = min(abs(wire_data(i,1:peak_idx(i))));
%             [~,t_1(i)] = min(abs(wire_data(i,peak_idx(i):valley_idx(i))));
%             accelerations(i,1) =  mean(diff(wire_data(i,1:peak_idx(i)),2),2);
%             accelerations(i,2) =  mean(diff(wire_data(i,peak_idx(i):end),2),2);
%         end