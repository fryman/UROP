%'\\chunky.mit.edu\v7_rodent\alexander friedman\clustered and unclastered data\rat1\2011-05-19_18-35-29_rat1_25 comb15\tt8pl6.ntt';
%'\\chunky.mit.edu\v7_rodent\alexander friedman\clustered and unclastered data\rat1\trains\2011-06-16_20-44-04_rat1_41_2combo15\tt11dms3.ntt';
%'\\chunky.mit.edu\v7_rodent\alexander friedman\clustered and unclastered data\rat1\2011-05-02_19-11-08_rat1_12_tr30';

%pathToFoo = 'F:\Aaron_Fryman\Data\';
foo = 'Data\rat1ses17tt23.ntt';
%fullFoo = fullfile(pathToFoo, foo);

[TimeStamps, DataPoints] = Nlx2MatSpike_411(foo, [1 0 0 0 1], 0, 1, []);
dp = reshape(DataPoints, 128, []);
% figure; plot(dp(:,1:100))
% size(dp)
save('rat1ses17tt23');