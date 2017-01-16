%% network expansion algorithm
% author: Joshua Goldford
% date: 11-6-2015

% inputs: R := reactant binary matrix (m x n)
%         P := product binary matrix (m x n)
%         x := starting point (binary vector (m x 1)
%         b := total number of reactant vector (n x 1) (can be computed via
%         row sums on R)
% output: out := output structure containing following fields:
%              x := binary vector representing presence/absence of metabolites 
%                   in converged network
%              y := binary vector representing presence/absence of reactions 
%                   in converged network
%              X := (metabolite accumulation matrix) binary mxr matrix 
%                   representing the presence/absence of metabolites at 
%                   iteration r
%              Y := (reaction accumulation matrix) binary mxr matrix 
%                   representing the presence/absence of reactions at 
%                   iteration r


function [out] = netExp(R,P,x,b)

% initialize variables:
% find the total number of metabolites in the seed set
k = sum(x);
% initialize previous iteration count of metabolites in the network
k0 = 0;
% iteration 1 consistes of the seed set
X = x;
% initialize reaction accumulation matrix
Y = [];

% while the metabolite set has not converged
while k > k0 
    % update previous number of metabolites
    k0 = k;
    % R'x ==> represnts the number of present metabolites in the 
    % network within each reaction; if this isequal to the total 
    % number of metabolites in each reaction, then the reaction 
    % is added to the network
    y  = double(R'*x == b);
    
    % P*y > 0 ==> represents the vector of reactions that produce 
    % metabolite i. (i in 1:m).  If this is >0, 
    % then that metabolite is producable 
    x_new = double(P*y > 0);
    % add to previous set of meatabolites (only needed to retain seed set)
    x = double(x | x_new);
    % find new total number of metabolites in network
    k = sum(x);
    
    % append accumulation matricies
    X = [X,x];
    Y = [Y,y];
end
    
% parse variables into output structure    
out.x = x;
out.y = y;
out.X = X;
out.Y = Y;
   


end

