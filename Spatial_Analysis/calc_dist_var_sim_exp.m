function varargout = calc_dist_var_sim_exp(INDist,KNDist,nsamp,tp,name)
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
INDist=struct2cell(INDist);
KNDist=struct2cell(KNDist);

% Inter-Nuclei distance
for i=1:nsamp

    INDist1 = INDist{i};
    INDist2 = INDist{i+nsamp};
    
    for j=1:length(tp)

    maxINDist = max(max(INDist1.tp{j}),max(INDist2.tp{j}));
    LstCount = maxINDist/100;

    pts = 0:LstCount:maxINDist;

    [f1,~] = ksdensity(INDist1.tp{j}, pts);
    [f2,~] = ksdensity(INDist2.tp{j}, pts);

    csIN = cosine_sim(f1,f2);
    cos_sim_IN = [cos_sim_IN;i,nsamp,csIN];
    
    end
end

% Nuclei Nearest Neighbor distance
for i = 1:nsamp

    NNDist1 = KNDist{i};
    NNDist2 = KNDist{i+nsamp};
    
    for j=1:length(tp)

    maxNNDist = max(max(NNDist1.tp{j}),max(NNDist2).tp{j});
    LstCount = maxNNDist/100;

    pts = 0:LstCount:maxNNDist;

    [f1,~] = ksdensity(NNDist1.tp{j}, pts);
    [f2,~] = ksdensity(NNDist2.tp{j}, pts);

    %[kl_dist,shannon_nn] = KLDiv(f1,f2);
    %KLdiv_NN = [KLdiv_NN;i,j,kl_dist,shannon_nn];
    csNN = cosine_sim(f1,f2);
    cos_sim_NN = [cos_sim_NN;i,i+nsamp,csNN];
    
    end   
end

writematrix(cos_sim_IN,[dirVar '/' 'Var_INDist_' name '.txt']);
writematrix(cos_sim_NN,[dirVar '/' 'Var_KNDist_' name '.txt']);