function [L2Err,L1Err,SE,omp_x,omp_supp,avg_FP,avg_FN] = OMP_algo(A,x,b,eps)
%Orthagonal Matching Persuit
%   Detailed explanation goes here
new_supp = [];
res = b;

iteration = 100;
for k = 1:iteration
    % Sweep Stage
    new_res = res;
%     err = new_res'*new_res - (A'*new_res).^2;
    proj = (A'*new_res);
    [~,max_pos] = max(abs(proj));
    
    % Update Support
    new_supp = union(new_supp,max_pos);
    
    % Update Provisonal Solution    
    x_s = pinv(A(:,new_supp))*b;
    new_x = zeros(size(x,1),1);
    new_x(new_supp,1) = x_s;
    
    % Update residual
    res = b - A*new_x;
    
    % Stoping rule
    if norm(res) < eps
        break;
    end
end
omp_supp = new_supp;
omp_x = new_x;

[L2Err,L1Err,SE,avg_FP,avg_FN] = matric_cal(new_x,x);


end

