clear; close all;
load('Data\rat1ses17tt8.mat');

N = 3;

newDP1 = dp(1:32,:)';
data_mean1 = mean(newDP1);
[data_covar1, sc1, la1] = princomp(newDP1);
peaks1 = max(newDP1, [], 2);
PCEL1 = sc1(:,1:N)';

newDP2 = dp(33:32*2,:)';
data_mean2 = mean(newDP2);
[data_covar2, sc2, la2] = princomp(newDP2);
peaks2 = max(newDP2, [], 2);
PCEL2 = sc2(:,1:N)';

newDP3 = dp(32*2+1:32*3,:)';
data_mean3 = mean(newDP3);
[data_covar3, sc3, la3] = princomp(newDP3);
peaks3 = max(newDP3, [], 2);
PCEL3 = sc3(:,1:N)';

newDP4 = dp(32*3+1:32*4,:)';
data_mean4 = mean(newDP4);
[data_covar4, sc4, la4] = princomp(newDP4);
peaks4 = max(newDP4, [], 2);
PCEL4 = sc4(:,1:N)';

peaks = [peaks1'; peaks2'; peaks3'; peaks4'];
allPC = [PCEL1; PCEL2; PCEL3; PCEL4];

inp = [TimeStamps; peaks; allPC];

[bestk, bestpp, bestmu, bestcov, dl, countf, indicators, numK] = fit_gauss_mixtures(inp, 1, 1, 10, 1e-10, 1e-4);
%save('em_clustered.mat','bestcov','bestk','bestmu','bestpp','exptsc');