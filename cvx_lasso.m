function [x]=cvx_lasso(b,A,lambda)
% this is the function using cvx to solve the LASSO problem
%Inputs - b.............measurements
%       - A.............measuring matrix
%       - lambda........parameter in the objective function
%Output - x.............reconstructed signal

n=size(A,2);
cvx_begin quiet
     variable x(n)
     minimize (sum(square_abs(b-(A*x)))+(2*lambda*norm(x,1)))
cvx_end