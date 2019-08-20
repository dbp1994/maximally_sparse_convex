function [L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(new_x,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% omp_supp = supp;
% omp_x = new_x;
eps_se = 1e-3;
s_xn = zeros(size(x,1),1);
s_x = zeros(size(x,1),1);
% Support Error
s_xn(find(abs(new_x)>eps_se)) = 1;
s_x(find(abs(x)>eps_se)) = 1;
SE = sum(abs(s_xn - s_x));

% x_idx = zeros(size(x,1),1);
% newx_idx = zeros(size(x,1),1);
% % FP = 0;
% % FN = 0;
% 
% newx_idx(new_supp) = 1;
% x_idx(supp_x) = 1;
% x_idx = [x_idx' zeros(1,length(x)-length(x_idx))];
% newx_idx = [newx_idx' zeros(1,length(x)-length(newx_idx))];

FP_FN = s_x - s_xn;
FP = find(FP_FN == -1);
if isempty(FP)
    avg_FP = 0;
else
    
    avg_FP = length(FP);
end
FN = find(FP_FN == 1);
if isempty(FN)
    avg_FN = 0;
else
    
    avg_FN = length(FN);
end


% L2Err = norm(normc(x) - normc(new_x));
L2Err = norm((x) - (new_x));
L1Err = norm(((x) - (new_x)),1);
% max_supp = max([length(new_supp),length(supp_x)]);
% distErr = (max_supp - length(intersect(supp_x,new_supp)))/max_supp;
end

