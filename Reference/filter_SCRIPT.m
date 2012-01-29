close all;

number_of_clusters = max(tt8pl6(:,2)) + 1;
clusters_indices = {1,number_of_clusters};
for i=1:number_of_clusters
    clusters_indices{1,i} = find(tt8pl6(:,2) == (i-1));
end

nrows = size(tt8pl6,1);
hw = hanning(17);
number_of_samples = 300;
numsamples = number_of_samples;

% noise_Nlenergy = zeros(nrows,4);
smoothwave_noise = zeros(4,numsamples);

for wire=1:4
    startIdx = 3+(wire - 1)*32 + 1;
    endIdx = 2+32*wire;
    for i=1:numsamples
        %         noise_Nlenergy(i,wire) = nonlinear_energies(clusters_indices{1,1}(i),wire);
        waveform = tt8pl6(clusters_indices{1,1}(i),startIdx:endIdx)...
            - mean((tt8pl6(clusters_indices{1,1}(i),startIdx:endIdx)));
        abswave = abs(waveform);
        smoothwave_noise(wire,i) = max(conv(abswave, hw, 'valid'));
    end
end

smoothwave_noise_max_wire = max(smoothwave_noise);

% cluster_idx=2;
for cluster_idx=2:number_of_clusters
    % cluster_NLenergy = zeros(nrows,wire);
    smoothwave_cluster = zeros(4,numsamples);
    
    for wire=1:4
        startIdx = 3+(wire - 1)*32 + 1;
        endIdx = 2+32*wire;
        for i=1:numsamples
            waveform = tt8pl6(clusters_indices{1,cluster_idx}(i),startIdx:endIdx)...
                - mean((tt8pl6(clusters_indices{1,cluster_idx}(i),startIdx:endIdx)));
            abswave = abs(waveform);
            smoothwave_cluster(wire,i) = max(conv(abswave, hw,'valid'));
        end
    end
    figure; hold on;

    smoothwave_cluster_max_wire = max(smoothwave_cluster);
    plot(1:numsamples,smoothwave_cluster_max_wire,'g')
    plot(1:numsamples,smoothwave_noise_max_wire, 'b')
end
%
%         plot(1:numsamples,noise_Nlenergy(1:numsamples,wire),'g');
%         plot(1:numsamples,cluster_NLenergy(1:numsamples,wire),'b');