function [x] = conjugate_grad(b,A,tol)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x = zeros(size(b,1),1);
r = b - A * x;
p = r;
rsold = r' * r;

while(1)
    Ap = A * p;
    alpha = rsold / (p' * Ap);
    x = x + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;
    if sqrt(rsnew) < tol
        break;
    end
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
end
end


