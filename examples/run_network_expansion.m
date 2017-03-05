% Author: Joshua Goldford
% Date: March 5th, 2017
% Title: Example script to run network expansion
% Description: This script will run the thermodynamically unconstrained
% network expansion algorithm presented for the following seed set: Water,
% Co2, Ammonia, Acetate, Formate, Hydrogen Sulfide, Bicarbonate, Nitrogen
% Gas.  The final variables

% Output are cell arrays called 'rxns' and 'mets' assinged to KEGG reaction
% and metabolites IDs, respectively, producable after network expansion.

clear
% PUT PATH TO networkExpansion DIRECTORY HERE
homedir = '/Users/joshuagoldford/Dropbox/phd/github/networkExpansion';
load([homedir,'/data/networks/network_filtered.mat']);

% add the network expansion function to your path
addpath([homedir,'/lib']);

% define the seed set
seed_set = { ...
    'C00001', ... % Water
    'C00011', ... % CO2
    'C00014', ... % Ammonia
    'C00033', ... % Acetate
    'C00058', ... % Formate
    'C00283', ... % Hydrogen sulfide
    'C00288', ... % Bicarbonate
    'C00697', ... % Nitrogen gas
    };


% block o2 production by defining a set of irreversible reations
% ensure all that produce o2 are in the o2 producing direction
o2  = find(strcmp('C00007',{network.compounds.cid}));

%find all reactions that produce o2
r = find(network.S(o2,:)' > 0);
% switch the stoichiometric coefficients to make these o2 consuming
% reactions
network.S(:,r) = -network.S(:,r);
% definte all reactions that consume o2 as irreversible
irrev = (network.S(o2,:)' < 0);


%% For each reversible reactions, create unidirectional forward and backward reactions, and create a new stoichiometric matrix
% convert stoichiometrix matrix into nx1 cell
S = num2cell(network.S,1);

% find (ir)reversible reactions
rev = find(~irrev);
ir = find(irrev);
rxns = network.rxns;

% for each reversible reaction, a forward and backward unidirectional
% reaction
for i =1:length(rev)
    % convert single stoichiometric vector into both forwar and backward
    % reactions
    S{rev(i)} = [S{rev(i)},  -S{rev(i)}];
    % update reaction names with _f for forward and _b for backward
    rxns{rev(i)} = {strcat(rxns{rev(i)},'_f'),strcat(rxns{rev(i)},'_b')};
   
end

% expand S-cell into S-matrix
S = cell2mat(S);
% expand reactions cell array
rxns = [rxns{:}]';

% Define input variables for network expansion.

% reactant matrix is 1 is s(i,j) < 0, 0 otherwise
input.R = double(S < 0);
% product matrix is 1 is s(i,j) > 0, 0 otherwise
input.P = double(S > 0);
% b(j) is the sum of all reatants for reation j
input.b = sum(input.R)';
% append the compound structure from the original network
input.compounds = network.compounds;
% append updated reaction network
input.rxns = rxns;


%% perform network expansion
% find the index of the seed_set
[~,seed] = intersect({input.compounds.cid},seed_set);

% map seet set to a vector 
x0  = zeros(length(input.compounds),1);
x0(seed) = 1;

% run network expansion
out = netExp(input.R,input.P,x0,input.b);

%% parse results 
% define helper anonymous functions for network expansion
getRxns = @(res) unique(cellfun(@(z) z(1:6), input.rxns(~~(res.y)),'uni',0));
getMets = @(res) {input.compounds(~~(res.x)).cid};
getMetStructs = @(res) input.compounds(~~(res.x));

% assign metabolites ids to variables
mets = getMets(out);
rxns = getRxns(out);
%metStructs = getMetStructs(out);

clearvars -except mets rxns 


