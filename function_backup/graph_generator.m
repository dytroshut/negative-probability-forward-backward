function [A,G] = graph_generator(n,x)

% n is the number of nodes.
% probability of edge between a pair of nodes

while 1
 
%Create a random adjacency matrix with n nodes.
R = exprnd(x,[n,n]);
R(R>1) = 1;
A = round(R);
A = triu(A) + triu(A,1)';
A = A - diag(diag(A));
%Since we have the adjacency matrix A, we can construct the graph G.
G1 = graph(A); %Get the incidence matrix I of the graph G

if length(bfsearch(G1,1)) == n
   G = G1;
   break;
end
   
end