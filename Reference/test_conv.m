%% Needs loaded: conv_filter output, process_data_alg output (both are Algorithm_Inputs).
close all;

axis1=2; axis2=3;

figure; hold on;
% color_values = smoothwave_max_wire/max(smoothwave_max_wire);
color_values = sigmoid_rating;
R = zeros(size(color_values,2),1);
G = 1 - color_values';
B = color_values';
cmap = [R G B];

scatter(peak_values(:,axis1), ...
    peak_values(:, axis2), ...
    0.5, cmap, '.');
whitebg
hold off;