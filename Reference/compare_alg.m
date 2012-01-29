function compare_alg( file_without_extension, axes )
%COMPARE_ALG Summary of this function goes here
%   Detailed explanation goes here
close all;
fwe = file_without_extension;

pathToFoo = 'C:\Users\Alexander\Desktop\Aaron Fryman\Results\Algorithm\rat1ses12tt8_exps'; %'F:\AaronF\'; %
foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
load(foo);

data = inp;

% data = [inp(1:5,:); -1*inp(6:end,:)];
[~, assignments] = max(optIndicators);
graph_clustered_data([assignments-1; data], optk, optIndicators, sigmoid_rating, axes);

% pathToFoo = 'F:\AaronF\Data\Manual\';
% foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
% loadedData = struct2cell(load(foo));
% 
% clusteredMat = loadedData{1}';
% 
% manualAssignments = clusteredMat;
% 
% k = max(manualAssignments(1,:), [], 2) + 1;
% graph_clustered_data(manualAssignments, k, ones(1, size(manualAssignments, 2)), axes);


end

