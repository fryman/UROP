function [ result_idx ] = ga_init_fn( data_points )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

num_components = 15;

[~,num_samples] = size(data_points);
rand_indices = randperm(num_samples);
result_idx = rand_indices(1:num_components);

end