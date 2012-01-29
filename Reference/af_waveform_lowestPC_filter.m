function [ cluster_ratings ] = af_waveform_lowestPC_filter( waveform_data, assignments, number_of_clusters, peak_vals )
%AF_WAVEFORM_FILTER Summary of this function goes here
%   waveform_data must not contain timestamps!

cluster_ratings = zeros(1,number_of_clusters);
min_cluster_size = sum(assignments==1);

for cluster=1:number_of_clusters
    min_cluster_size = min(min_cluster_size, sum(assignments==cluster));
end

for cluster=1:number_of_clusters
%     cluster_data = waveform_data(assignments==cluster, :);
    cluster_data = peak_vals(assignments==cluster,:);
    rand_indices = randperm(min_cluster_size);
    cluster_variances = var(cluster_data(rand_indices,:))/mean(cluster_data(rand_indices,:));
%     cluster_variances = var(cluster_data);
    cluster_ratings(cluster) = mean(cluster_variances);
end

end