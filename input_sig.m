function [r] = random_sig(interval,N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
r = interval(1) + (interval(2)-interval(1))*rand(N,1);

end

