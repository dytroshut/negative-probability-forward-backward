function x = FBsolver(a,b,c)

% a,b are force to be the column vector with same size
% x = log( ( -c+sqrt(c.^2 + 4.*a.*b) )./ (2.*b) );
x = log( ( -c+sqrt(c.^2 + 4.*a.*b) )./ (2.*b) );
% x(x==-Inf) = Inf;
% x(x==Inf) = -Inf;

end