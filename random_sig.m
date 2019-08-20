function [r,K,supp] = random_sig(intrvl_t,intrvl_x,K,N)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

t = randi([intrvl_t(1),intrvl_t(2)],K,1);
idx = 0;
r = zeros(N,1);
for i = 1:length(t)
    idx = idx + t(i);
    if idx >= N
        break;
    end
    r(idx) = intrvl_x(1) + (intrvl_x(2)-intrvl_x(1))*rand();
    supp(i) = idx;
end
K = length(supp);
end

