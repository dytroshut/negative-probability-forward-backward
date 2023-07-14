function [A] = negative_generator(n,x)

R = exprnd(x,[n,n]);
R(R>1) = 1;
A = round(R);
A = triu(A) + triu(A,1)';
A = A - diag(diag(A));

end