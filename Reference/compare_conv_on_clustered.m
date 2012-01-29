%% Needs loaded: conv_filter output, manual assignment matrix.

% fwe = 'conv_1-12-9';
% pathToFoo = 'C:\Users\Alexander\Desktop\Aaron Fryman\Data\Algorithm_Inputs';
% foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
% load(foo);

fwe2 = '1-12-9';
pathToFoo2 = 'C:\Users\Alexander\Desktop\Aaron Fryman\Results\Manual';
foo2 = fullfile(pathToFoo2, strcat(fwe2, '.mat'));
load(foo2);

close all;

number_of_clusters = max(manual_output(:,2)) + 1;
clusters_indices = {1,number_of_clusters};

for i=1:number_of_clusters
    clusters_indices{1,i} = find(manual_output(:,2) == (i-1));
end

nrows = size(manual_output,1);
number_of_samples = 500;
num = number_of_samples;

smoothwave_noise = smooth_max(clusters_indices{1,1});

for cluster_idx=2:number_of_clusters
    smoothwave_cluster = smooth_max(clusters_indices{1,cluster_idx});
    figure; hold on;

    plot(1:num,smoothwave_cluster(1:num),'g');
    plot(1:num,smoothwave_noise(1:num), 'b');
    title(strcat('CLUSTER ',num2str(cluster_idx-1)));
    hold off;
end