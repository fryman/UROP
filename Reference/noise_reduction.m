%NO LONGER USED
%I think this was going to be some form of overclustering using kmeans.
%Failed experiment.

%Parameters
numPointsPerComp = 500;

fwe = 'rat1ses12tt8';
pathToFoo = 'F:\AaronF\Data\Algorithm\'; %'C:\Users\Alexander\Desktop\Aaron Fryman\Data'; %
foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
load(foo);

N = 3;
newDP = dp';
nrows = size(newDP, 1);
score = zeros(nrows, 128);
peaks = zeros(nrows, 4);
prin_comps = zeros(12, nrows);
bad_wires = [0 0 0 0];

for wire = 1:4
    startIdx = (wire - 1)*32 + 1;
    endIdx = 32*wire;
    if mean(newDP(:, startIdx:endIdx)) ~= zeros(1, 32)
        [~, score(:, startIdx:endIdx)] = princomp(newDP(:, startIdx:endIdx));
        peaks(:, wire) = max(newDP(:, startIdx:endIdx), [], 2);
        prin_comps((wire - 1)*3 + 1:wire*3, :) = score(:, startIdx:startIdx + N - 1)';
    else
        bad_wires(wire) = 1;
    end
end

if any(bad_wires)
    for badIdx = find(bad_wires)
        startIdx = (badIdx - 1)*3 + 1;
        endIdx = 3*badIdx;
        peaks(:, badIdx) = [];
        prin_comps(startIdx:endIdx, :) = [];
    end
end

inp = [TimeStamps; peaks'; prin_comps];

initK = floor(size(inp, 2)/numPointsPerComp);
IDX = kmeans(inp(2:3,:)', initK); %Which dimensions to kmeans cluster?