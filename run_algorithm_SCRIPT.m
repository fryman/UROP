% close all; clear;
% load('C:\Users\Alexander\Desktop\Aaron Fryman\Data\Algorithm_Inputs\rat1ses12tt8.mat')
% load('C:\Users\Alexander\Desktop\Aaron Fryman\Data\Algorithm_Inputs\conv_1-12-8.mat')
% 
% %PARAMETERS
% baseline_noise_threshold = 275;
% data = [TimeStamps; peak_values'; prin_comps];
% num_runs = 1;
% weight_threshold_for_mean_randomization_indices = 0.5;
% 
% noise_idx = (smoothwave < baseline_noise_threshold);
% noise_feature_data = data(:,noise_idx);
% cluster_feature_data = data(:,noise_idx==0);
% input_weights = sigmoid_rating(noise_idx==0);
% 
% 
% while True
% %     optk_results = zeros(1,num_runs);
% %     indicators_results = cell(1,num_runs);
% %     normalized_indicators_results = cell(1,num_runs);
%     for i=1:num_runs
%         disp('RUN NUMBER'), disp(i)
%         [number_of_components, mixture_probabilities, means, covariances, description_length, ~, indicators, ~, normalized_indicators] = fit_gauss_mixtures(cluster_feature_data, input_weights, 1, 15, 1e-10, 1e-4, weight_threshold_for_mean_randomization_indices);
% 
% %         optk_results(i) = number_of_components;
% %         indicators_results{i} = optIndicators;
% %         normalized_indicators_results{i} = bestNormIndicators;
%         [~, assignments] = max(indicators);
%         reconstructed_inputs = [assignments zeros(1,size(noise_feature_data,2)); cluster_feature_data noise_feature_data]';
%         reconstructed_noise_rating = [input_weights sigmoid_rating(noise_idx == 0)]';
%         max_filter_values_per_cluster = af_waveform_filter(dp', assignments, number_of_components);
%         LRatios = af_rateClusters([reconstructed_inputs(:,1) reconstructed_inputs]);
%     end
% end
% 
% %%
% for idx=1:size(LRatio_clusters,2)
%     cluster_data = reconstructed_input(full_assignments==LRatio_clusters(idx),:);
%     LRatio_data((1 + offset):(size(cluster_data,1) + offset), 1:end) = cluster_data;
% end
% %%
% offset = 0;
% for idx=1:size(LRatio_clusters,2)
%     cluster_data = reconstructed_input(full_assignments==LRatio_clusters(idx),:);
%     if idx > 1
%         offset = offset + sum(full_assignments==LRatio_clusters(idx-1)); %we want to insert at the end of the previous data. offset tracks cursor.
%     end
%     LRatio_data((1 + offset):(size(cluster_data,1) + offset), 3:end) = cluster_data;
% end
%         
% figure; plot(description_length,'LineWidth',2); title('Description Length');
 

% graph_clustered_data(reconstructed_data, optk, bestNormIndicators, reconstructed_noise_rating, [3 4]);
% graph_clustered_data([assignments-1; inp(:,noise)], optk, bestNormIndicators, sigmoid_rating(noise), [3 4]);
%%OVERCLUSTERING

%IMPORTANT!!: this is the code that removes the baseline noise cluster!
%This is will need to be made into a more general, useable function, or
%appended to an already existing one.
% clear;load('input.mat')
baseline_noise_threshold = 0.1;
weight_threshold_for_mean_randomization_indices = 0.5;
waveform_filter_threshold = 21;

baseline_noise_idx = (conv_filter_rating < baseline_noise_threshold);
noise_feature_data = input(:,baseline_noise_idx);
cluster_feature_data = input(:,baseline_noise_idx==0);
loop = 1;

input_weights = conv_filter_rating(baseline_noise_idx==0);
reconstructed_noise_rating = [input_weights conv_filter_rating(baseline_noise_idx)];
reconstructed_input = [cluster_feature_data noise_feature_data];
reconstructed_waveform_data = [waveform_data(baseline_noise_idx==0, :); waveform_data(baseline_noise_idx, :)];
reconstructed_noise_feature_data = noise_feature_data;
reconstructed_cluster_feature_data = cluster_feature_data;

clear baseline_noise_idx noise_feature_data cluster_feature_data

%% overclustering
[number_of_components, ~, ~, ~, ~, ~, indicators, ~, normalized_indicators] = ...
    fit_gauss_mixtures(reconstructed_cluster_feature_data, input_weights, 30, 50, 1e-10, 1e-4, weight_threshold_for_mean_randomization_indices); %OVERCLUSTERING
[~, assignments] = max(indicators);
full_assignments = [assignments zeros(1,size(reconstructed_noise_feature_data,2))];


max_filter_values_per_cluster = af_waveform_filter(reconstructed_waveform_data, full_assignments, number_of_components);
disp(number_of_components)



noise_idx = (full_assignments==0);
noise_cluster_idx = find(max_filter_values_per_cluster<21);
for i=1:size(noise_cluster_idx,2)
    noise_idx = logical(noise_idx + (full_assignments==noise_cluster_idx(i)));
end

input_weights = reconstructed_noise_rating(noise_idx==0);
reconstructed_noise_rating = [input_weights reconstructed_noise_rating(noise_idx)];

reconstructed_cluster_feature_data = reconstructed_input(:, noise_idx==0);
reconstructed_noise_feature_data = reconstructed_input(:, noise_idx);
reconstructed_waveform_data = [reconstructed_waveform_data(noise_idx==0, :); reconstructed_waveform_data(noise_idx, :)];
reconstructed_input = [reconstructed_cluster_feature_data reconstructed_noise_feature_data];

%normal clustering
[number_of_components, ~, ~, ~, ~, ~, indicators, ~, normalized_indicators] = ...
    fit_gauss_mixtures(reconstructed_cluster_feature_data, input_weights, 1, 15, 1e-10, 1e-4, weight_threshold_for_mean_randomization_indices); %OVERCLUSTERING
[~, assignments] = max(indicators);
full_assignments = [assignments zeros(1,size(reconstructed_noise_feature_data,2))];
    
graph_clustered_data(reconstructed_input, full_assignments, number_of_components, normalized_indicators, reconstructed_noise_rating', [3 4]);
