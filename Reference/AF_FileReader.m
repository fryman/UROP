p = path;
path(p,'F:\dgibson\nlx_411\');
MyFolder='D:\Tradeof_project_Data_in_Work\Rat1\2011-05-07_17-40-46_rat1_17_cost_cost3_is_4min';
muFile='tt3pl1.ntt';
myfile=fullfile(MyFolder, muFile);
%myfile = '\\chunky.mit.edu\v7_rodent\alexander friedman\clustered and unclastered data\rat1\unclasterd\2011-05-19_18-35-29_rat1_25 comb15\tt8pl6.ntt';
DataPoints = Nlx2MatSpike_411(myfile, [0 0 0 0 1], 0, 1, []);
dp = reshape(DataPoints, 128, []);
figure; plot(dp(:,1:100))
size(dp)
%save('2011-05-19_18-35-29_rat1_25 comb15-tt8pl6', 'dp');

path(p);