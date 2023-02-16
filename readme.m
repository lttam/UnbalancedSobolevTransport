% This is the code for the paper
% Title: "Scalable Unbalanced Sobolev Transport for Measures on a Graph"
% Authors: Tam Le, Truyen Nguyen, Kenji Fukumizu
% Published at AISTATS 2023



% ***** Data: e.g., 'twitter.mat' for the TWITTER dataset
% Link: https://www.dropbox.com/s/nhoor4jnvfd0xlk/twitter.mat?dl=0



% ********** Third-party toolboxes

% *** -- For building graph (G_Log / G_Sqrt) from support data points
% + we use the third-party toolbox
% Link: https://github.com/lttam/SobolevTransport

% *** -- For Sinkhorn-based UOT (and the given graph computed by the above
% third-party toolbox)
% we use the third-party toolbox
% Link: https://github.com/gpeyre/2017-MCOM-unbalanced-ot

% *** -- Note: the dAlpha of EPT on a tree is a special case of our
% proposed Unbalanced Sobolev Transport (UST) (see Proposition 5.3 i)).
% Therefore, we use our implementation of UST to compute dAlpha (of EPT on
% a tree) when the graph G is a tree.



% ********** Given the graph (computed by the third-party toolbox)

% *** -- Compute distance matrices for Unbalanced Sobolev transport and EPT on trees
% + compute_UnbalancedSobolevTransport_vUS: compute the distance matrix for unbalanced Sobolev
% transport (and some variants)
%
% + compute_EPT_dAlpha_randTree: compute the distance matrix for EPT on
% trees randomly sampled from the given graph (and some variants of EPT)

% *** -- Note:
% The code uses Graph and Network Algorithms toolbox from MATLAB. (e.g., Dijkstra
% algorithm for shortest path from a source point to a destination set of
% points.)


