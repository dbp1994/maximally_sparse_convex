function [A,B] = gen_toeplitz(a,b,N)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
col_A = [a zeros(1,N-length(a))];
col_B = [b zeros(1,N-length(b))];

row_A = [a(1) zeros(1,N-1)];
row_B = [b(1) zeros(1,N-1)];

A = toeplitz(col_A,row_A);
B = toeplitz(col_B,row_B);

end

