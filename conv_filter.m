function [ sigmoid_rating, smoothwave_max_wire ] = conv_filter( file_without_extension )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

fwe = file_without_extension;
pathToFoo = 'C:\Users\Alexander\Desktop\Aaron Fryman\Data\';
foo = fullfile(pathToFoo, strcat(fwe, '.mat'));
load(foo);

hw = hanning(17);
[~, ncols] = size(dp);
smoothwave = zeros(4,ncols);

for wire=1:4
    startIdx = (wire - 1)*32 + 1;
    endIdx = 32*wire;
    for i=1:ncols
        waveform = dp(startIdx:endIdx,i)...
            - mean((dp(startIdx:endIdx,i)));
        abswave = abs(waveform);
        smoothwave(wire,i) = max(conv(abswave, hw, 'valid'));
    end
end

smoothwave_max_wire = max(smoothwave);
sigmoid_rating = 1./(1+exp(5*10^-2*(-1*smoothwave_max_wire+300)));

end