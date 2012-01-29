load('noisebyenergy.mat')
% global individual_id;
% individual_id = 0;

popsize = 6;
gennum = 100;
num_elite = 1;
crossover = .8;

num_variables = 15;
num_samples = length(reconstructed_cluster_feature_data);
randindex = randperm(num_samples);
initPop = randindex(1:num_variables);

gaoptions = gaoptimset('PopulationSize',popsize,'Generations',gennum,'EliteCount',num_elite,...
    'PlotFcns', @gaplotbestf, 'InitialPopulation', initPop, 'CrossoverFraction',crossover,...
    'UseParallel','always');
lb=ones(1,num_variables);
ub=repmat(size(reconstructed_cluster_feature_data,2),[1 num_variables]);
IntCon = 1:num_variables;

if matlabpool('size') == 0
    matlabpool
end

% for x=1:5
parfor x=1:8;
    FitnessFcn = @(rand_indices)af_run_algorithm(rand_indices, reconstructed_input',...
        reconstructed_cluster_feature_data, input_weights, reconstructed_noise_feature_data,...
        reconstructed_waveform_data, reconstructed_noise_rating, 'sum', x);
    
    [xGA, fGA, exitflag, output] = ga(FitnessFcn, num_variables,[],[],[],[],lb,ub,[],IntCon,gaoptions);
    results = struct('BestIndividual', xGA, 'BestFitnessVal',fGA,'ExitFlag',exitflag,'Output',output);
    foo = sprintf('wf22_sum_%d.mat', x); %fullfile('/home/fryman/Results/Genetic_Algorithm/wff22/', sprintf('wf22_sum_%d.mat', x));
    parsave(foo, 'gaoptions', 'results');
end
matlabpool close;