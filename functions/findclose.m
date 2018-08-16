function [ closest_param ] = findclose( param, x, target )
%FINDCLOSE Summary of this function goes here
%   Detailed explanation goes here
index_train = 1:length(x);
[~,index] = max(index_train((x-target)<0));
closest_param = param(index);
end

