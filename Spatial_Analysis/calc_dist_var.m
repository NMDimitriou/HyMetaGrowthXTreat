function varargout = calc_dist_var(INDist,KNDist,name)
%% Computes variation between distributions using the Cosine similarity measure
% Input: 
%   INDist: Inter-Nuclei distances
%   KNDist: Nearest-Neighbour distances
%   name1, name2: the names for the output files
% Output: the files with the cosine similarity of the distributions
% Author: Nikolaos M. Dimitriou, 
% McGill University, 2020

dirVar='Distances/Variability'; % directory to save the variability between distances


if ~exist(dirVar, 'dir')
    mkdir(dirVar)
end


cos_sim_IN = [];
cos_sim_NN = [];
%INDist=struct2cell(INDist);
KNDist=struct2cell(KNDist);
%{
% Inter-Nuclei distance
for i = 1:length(INDist)-1
    for j = i+1:length(INDist)

    INDist1 = INDist{i};
    INDist2 = INDist{j};

    maxINDist = max(max(INDist1),max(INDist2));
    LstCount = maxINDist/100;

    pts = 0:LstCount:maxINDist;

    [f1,~] = ksdensity(INDist1, pts);
    [f2,~] = ksdensity(INDist2, pts);

    csIN = cosine_sim(f1,f2);
    cos_sim_IN = [cos_sim_IN;i,j,csIN];

    end
end
%}
% Nuclei Nearest Neighbor distance
for i = 1:length(KNDist)-1
    for j = i+1:length(KNDist)

    NNDist1 = KNDist{i};
    NNDist2 = KNDist{j};

    maxNNDist = max(max(NNDist1),max(NNDist2));
    LstCount = maxNNDist/100;

    pts = 0:LstCount:maxNNDist;

    [f1,~] = ksdensity(NNDist1, pts);
    [f2,~] = ksdensity(NNDist2, pts);

    %[kl_dist,shannon_nn] = KLDiv(f1,f2);
    %KLdiv_NN = [KLdiv_NN;i,j,kl_dist,shannon_nn];
    csNN = cosine_sim(f1,f2);
    cos_sim_NN = [cos_sim_NN;i,j,csNN];
    

    end
end

writematrix(cos_sim_IN,[dirVar '/' 'Var_INDist_' name '.txt']);
writematrix(cos_sim_NN,[dirVar '/' 'Var_KNDist_' name '.txt']);
