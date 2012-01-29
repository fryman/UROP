clear;clc;close all;

%% Creates 4 figures with random qi weights
% pcd = load_pcd('cv1.pcd');
% for clrs = 0:9
%     for it = 1:4
%         figure(it);
%         clf;
%         [K newPcd] = plot_weighted(pcd,clrs);
%     end
%     fprintf('cluster %d',clrs);
%     input(':');
% end
%% EM Algorithm
pcd = load_pcd('~/ar12/exp/model_cvs/shredded_wheat_cv1.pcd');
[K newPcd] = plot_weighted(pcd,0);
% for k = 1:numel(K)
%     if K(k) < .1
%         newPcd.data(k,:) = zeros(1,size(pcd.data,2));
%     else
%         K(k) = ceil(K(k)*10);
%         % fprintf('%d \n',K(k));
%         newPcd.data = [newPcd.data; repmat(newPcd.data(k,:),[K(k) 1])];
%     end
% end
% 
newPcd = populate_pcd_fields(newPcd.columns, newPcd.data);
y = [newPcd.X'; newPcd.Y'; newPcd.Z'];
%display(sum(K));
% K = ones(1000,1);
% K = [ones(1,250)*.00001, ones(1,250), ones(1,250), ones(1,250)*.00001];
% Y1 = [1:500,1:500];
% X1 = [100*ones(1,500), (-100)*ones(1,500)];
% figure(1);
% scatter(X1, Y1);
% X2 = normrnd(X1,10);
% figure(2);
% scatter(X2, Y1);
% y = [X2; Y1];

%[bestk, bestpp, bestmu, bestcov, dl, countf] = nd_gaus_mixtures_nonweighted(y,1,20,1e-15,1e-4,0);
[bestk, bestpp, bestmu, bestcov, dl, countf] = nd_gaus_mixtures_clean(y,K',1,20,1e-15,1e-4);
% [bestk, bestpp, bestmu, bestcov, dl, countf] = nd_gaus_mixtures(y,K',1,20,1e-15,1e-4,0);
figure
plot3(y(1,:),y(2,:),y(3,:),'o'); box on; grid on; hold on
% plot(y(1,:),y(2,:),'o'); box on; grid on; hold on
npoints = 300;
x = genmix(npoints,bestmu,bestcov,bestpp(1:end-1));
plot3(x(1,:),x(2,:),x(3,:),'r+');
% plot(x(1,:),x(2,:),'r+');
legend('original data','sampled from fit')
axis vis3d
axis equal

%cd F:/Aaron_Fryman
%load('Rat1Ses17tt3Cost.mat')
%[bestk, bestpp, bestmu, bestcov, dl, countf] = nd_gaus_mixtures_clean(dp,0,1,12,1e-15,1e-4);

% options = statset('Display','final');
% mat = [newPcd.X newPcd.Y newPcd.Z];
% gm = gmdistribution.fit(mat,3,'Options',options);
% idx = cluster(gm,mat);
% cluster1 = (idx == 1);
% cluster2 = (idx == 2);
% cluster3 = (idx == 3);
% %newPcd.clrIdx = idx;
% newPcd.R = double(cluster1);
% newPcd.B = zeros(length(pcd.X),1);%cluster2;
% newPcd.G = double(cluster3);
% figure(2);
% plot_pcd_color(newPcd);