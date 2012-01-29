function [config_rating] = af_run_algorithm( rand_indices, full_feature_data, cluster_feature_data, input_weights, noise_feature_data, wfdata, noise_rating, rating_mode, coeff_idx )
%AF_RUN_ALGORITHM Handles the FIT_GAUSS_MIXTURES fn and rates the
%cluster configuration
%   Detailed explanation goes here

% global base;
% global coeffs;
% global coeff_idx;
% coeff_idx = coeff_idx + 1;

weight_threshold_for_mean_randomization_indices = 0.5;
waveform_filter_threshold = 22;

[number_of_components, ~, ~, ~, ~, ~, indicators, ~, normalized_indicators] = ...
    fit_gauss_mixtures(cluster_feature_data(2:end,:), input_weights, 1, 15, 1e-10, 1e-4, weight_threshold_for_mean_randomization_indices, rand_indices);
if number_of_components > 1
    [~, assignments] = max(indicators);
else
    assignments = ones(size(indicators));
end
full_assignments = [assignments zeros(1, length(noise_feature_data))];

display(number_of_components);

% graph_clustered_data(full_feature_data', full_assignments, number_of_components, normalized_indicators, noise_rating', [3 4]);

LRatios = af_rateClusters([zeros(size(full_assignments')) full_assignments' full_feature_data]);
nonzero_LRatios = LRatios(LRatios ~= 0);
waveform_filter = af_waveform_filter(wfdata, full_assignments, number_of_components);
nonzero_wf_filter = waveform_filter(LRatios ~= 0);

cluster_sizes = zeros(1, number_of_components);
for cluster_idx=1:number_of_components
    cluster_sizes(cluster_idx) = sum(full_assignments == cluster_idx);
end

if sum(cluster_sizes < 100) == 0
    switch rating_mode
        case 'means' % this is the current standard
            config_rating = mean(LRatios);
            
        case 'sqrt' % sqrt of LRatio
            config_rating = mean(sqrt(nonzero_LRatios));
            
        case 'sqrt size limited' % sqrt of LRatio times avg size of noise cluster
            config_rating = mean(sqrt(nonzero_LRatios))*mean(cluster_sizes(nonzero_wf_filter < waveform_filter_threshold));
            
        case 'product'
            config_rating = mean(sqrt(nonzero_LRatios))*mean(LRatios)*mean(cluster_sizes(nonzero_wf_filter < waveform_filter_threshold));
            
        case 'sum'
            coeffs = [1 2 4 8 16 32 64 128];
            term1 = mean(sqrt(nonzero_LRatios))*mean(cluster_sizes(nonzero_wf_filter < waveform_filter_threshold))/10000;
            term2 = mean(LRatios);
            config_rating = coeffs(coeff_idx)*term1 + term2;
            
        case 'final'
            config_rating = mean(sqrt(nonzero_LRatios))*mean(cluster_sizes(nonzero_wf_filter < waveform_filter_threshold))/sum(nonzero_wf_filter >= waveform_filter_threshold);
            
        case 'sqrt over good' % sum of the sqrt of LRatio over the number of good clusters
            config_rating = sum(sqrt(nonzero_LRatios))/sum(nonzero_wf_filter >= waveform_filter_threshold);
            
        case 'no noisy' % only use good clusters as ratings
            config_rating = mean(nonzero_LRatios(waveform_filter >= waveform_filter_threshold));
            
        case 'weighted'
            best_clusters = (waveform_filter >= 30);
            worst_clusters = (waveform_filter < 21);
            weighted = ones(size(waveform_filter))*2 - best_clusters - worst_clusters*2;
            weighted_LRatios = weighted.*LRatios;
            config_rating = mean(weighted_LRatios(weighted_LRatios ~= 0));
            
        case 'sigmoid'
            % filter LRatio values through a single sigmoid that creates an
            % effective lower bound
            
        case 'filtered'
            high_wff_thresh = 30;
            mid_wff_thresh = 20;
            low_wff_thresh = 10;
            high_LR_thresh = 0.2;
            mid_LR_thresh = 0.4;
            low_LR_thresh = 0.6;
            
            artifact_wff_val = min(waveform_filter.*LRatios);
            if artifact_wff_val < 0.01
                waveform_filter = waveform_filter(waveform_filter ~= artifact_wff_val);
            end
            
            classification = struct('High', [], 'Middle', [], 'Low', []);
            classification.High = find(waveform_filter > high_wff_thresh);
            classification.Middle = find(waveform_filter > mid_wff_thresh - waveform_filter > high_wff_thresh);
            classification.Low = find(waveform_filter <= mid_wff_thresh);
            
            c = high_LR_thresh - high_LR_thresh/2;
            high_rating = mean(sigmf(LRatios(classification.High), [-10/c c]));
            c = (high_LR_thresh-mid_LR_thresh)/2 + mid_LR_thresh;
            mid_rating = mean(sigmf(LRatios(classification.Middle), [-10/c c]));
            c = (mid_LR_thresh-low_LR_thresh)/2 + low_LR_thresh;
            low_rating = mean(sigmf(LRatios(classification.Low), [-10/c c]));
            rating_terms = [high_rating, mid_rating, low_rating];
            config_rating = mean(rating_terms(logical(1-isnan(rating_terms))));
            
        case 'log'
            switch base
                case 'e'
                    config_rating = mean(log(nonzero_LRatios));
                case 10
                    config_rating = mean(log10(nonzero_LRatios));
                case 7.5
                    config_rating = mean(log10(nonzero_LRatios)/log10(7.5));
                case 5
                    config_rating = mean(log10(nonzero_LRatios)/log10(5));
                case 2
                    config_rating = mean(log2(nonzero_LRatios));
            end
    end
else
    config_rating = Inf;
end

end

        % THE FOLLOWING block of code can be used to seperate good and bad
        % cluster data
%         good_clusters = [0 find(max_filter_values_per_cluster >= waveform_filter_threshold)];
%         num_good_clusters = length(good_clusters);
%         good_clusters_length = 0;
%         for idx=1:num_good_clusters
%             good_clusters_length = good_clusters_length + sum(full_assignments==good_clusters(idx));
%         end
%         good_cluster_data = zeros(good_clusters_length, min(size(full_feature_data)) + 2);
%         offset = 0;
%         for idx=1:num_good_clusters
%             cluster_data = full_feature_data(full_assignments==good_clusters(idx),:);
%             num_cluster_spikes = size(cluster_data,1);
%             good_cluster_data((1 + offset):(num_cluster_spikes + offset), 2:end) = [repmat(idx-1, [num_cluster_spikes, 1]) cluster_data];
%             offset = offset + sum(full_assignments==good_clusters(idx)); %we want to insert at the end of the previous data. offset tracks cursor.
%         end
%         
%         bad_clusters = [0 find(max_filter_values_per_cluster < waveform_filter_threshold)];
%         num_bad_clusters = length(bad_clusters);
%         
%         bad_clusters_length = 0;
%         for idx=1:num_bad_clusters
%             bad_clusters_length = bad_clusters_length + sum(full_assignments==bad_clusters(idx));
%         end
%         
%         bad_cluster_data = zeros(bad_clusters_length, min(size(full_feature_data)) + 2);
%         offset = 0;
%         for idx=1:num_bad_clusters
%             cluster_data = full_feature_data(full_assignments==bad_clusters(idx),:);
%             num_cluster_spikes = size(cluster_data,1);
%             bad_cluster_data((1 + offset):(num_cluster_spikes + offset), 2:end) = [repmat(idx-1, [num_cluster_spikes, 1]) cluster_data];
%             offset = offset + sum(full_assignments==bad_clusters(idx)); %we want to insert at the end of the previous data. offset tracks cursor.
%         end