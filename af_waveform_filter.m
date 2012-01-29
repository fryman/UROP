function [ max_filter_values_per_cluster ] = af_waveform_filter( waveform_data, assignments, number_of_clusters )
%AF_WAVEFORM_FILTER Summary of this function goes here
%   waveform_data must not contain timestamps!

max_filter_values_per_cluster = zeros(1,number_of_clusters);

for cluster=1:number_of_clusters
    cluster_data = waveform_data(assignments==cluster, :);
    if size(cluster_data,1) > 1
        max_filter_values_per_cluster(cluster) = max(diff(std(cluster_data) - mean(cluster_data)));
    elseif size(cluster_data,1) == 0
        max_filter_values_per_cluster(cluster) = 0;
    else
        max_filter_values_per_cluster(cluster) = std(cluster_data) - mean(cluster_data);
    end
end

end