function graph_clustered_data( data, assignments, k, ~, conv_weights, axes )
%graph_clustered_data Summary of this function goes here
%   Detailed explanation goes here

axisLabels = {'TimeStamp', 'PK EL1', 'PK EL2', 'PK EL3', 'PK EL4', ...
    'PC1 EL1', 'PC2 EL1', 'PC3 EL1', 'PC1 EL2', 'PC2 EL2', 'PC3 EL2', ...
    'PC1 EL3', 'PC2 EL3', 'PC3 EL3', 'PC1 EL4', 'PC2 EL4', 'PC3 EL4', '1','2','3','4'};

axis1 = axes(1);
axis2 = axes(2);

for primaryComp=0:k
    figure
    hold on
%     whitebg
    
    secondaryData = data(:, (assignments ~= primaryComp));
    num_points = min(10000, size(secondaryData,2));
    
    randindex = randperm(size(secondaryData,2));
    randindex = randindex(1:num_points);
    scatter(secondaryData(axis1, randindex), ...
        secondaryData(axis2, randindex), ...
        0.01, [.5 .5 .5], '.');
    
%     color_values = normalized_indicators(assignments == primaryComp);
    color_values = conv_weights(assignments == primaryComp);
    R = zeros(size(color_values,1),1);
    G = 1 - color_values;
    B = color_values;
    cmap = [R G B];
%     cmap = [0.0 0.5 1.0]; %%%DEFAULT
    primaryData = data(:, (assignments == primaryComp));
    scatter(primaryData(axis1, :), ...
        primaryData(axis2, :), ...
        5, cmap, '.');
        
        
    set(gca,'FontName','Times','FontSize',22,'XLim',[-2*10^4, 4*10^4],'YLim',[0, 4*10^4]);
    xlabel(axisLabels(axis1));
    ylabel(axisLabels(axis2));
    placex = get(gca,'Xlim'); placey = get(gca,'Ylim');
    title(strcat('CLUSTER ',num2str(primaryComp)));
    text(placex(1)+0.2,placey(2)-0.2,sprintf('k=%d',k),...
        'FontName','Times','FontSize',22);
        
    drawnow
    hold off
    
%     figure
%     hist(color_values(1:end))
end

end

