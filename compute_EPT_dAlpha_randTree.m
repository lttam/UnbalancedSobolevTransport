% compute the distance matrix for EPT on
% trees randomly sampled from the given graph (and some variants of EPT)
% ***********************************************************************************

clear all
clc

maxKC = 100;
typeGGArray = {'RandLLE', 'RandSLE'};
ppArray = 1;

for typeGGID = 1:length(typeGGArray)  
        typeGG = typeGGArray{typeGGID};
        
% typeGG = 'RandLLE'; % log-linear #edges
% typeGG = 'RandSLE'; % sqrt-linear #edges

dsName = 'twitter';
nSS = 20; % #tree (average for Sobolev)

for ppID = 1:length(ppArray)
    
    pp = ppArray(ppID);

    disp(['-------- Type: ' typeGG ' && pp = ' num2str(pp) ' ---------']);
    
% DD_SS1, 5, 10, 20
load([dsName '_' num2str(maxKC) '_' typeGG '_Graph.mat']);

% nGG: number of vertices
randSArray = randperm(nGG);

wwGG = GG.Edges.Weight;

DD_SS = cell(nSS, 1); % Unbalanced Sobolev
varDD_SS = cell(nSS, 1); % variant of Unbalanced Sobolev

runTime_TreeSampling = zeros(nSS, 1);
runTime_Prep = zeros(nSS, 1);
runTime_Dist = zeros(nSS, 1);

%================
% OPT 
opt.a = 1;
opt.b = 1;
opt.lambda = 1;
opt.alpha = 0;
opt.c = 1;
opt.x_0 = 'r'; % root --> w(r) = a; [w(x) = d(r, x) + a]

for idSS = 1:nSS

    
% %     % ------- FOR EACH S0 (randomly choose) ---------
    tic
    disp('... randomly sampling tree');
    [minTR_EdgeID, minTR_EdgeWW] = RandomlySamplingTree(nGG, GG.Edges.EndNodes, GG.Edges.Weight);

    % minTR_EdgeID: nx2
    % minTR_EdgeWW: nx1
    
    minTRGG = graph(minTR_EdgeID(:, 1), minTR_EdgeID(:, 2), minTR_EdgeWW);

    s0=1;
    runTime_TreeSampling(idSS) = toc;

    tic
    disp(['...[' num2str(idSS) '] compute the tree path']);
    % tree path!!!
    [trPP, trDD, trEP] = shortestpathtree(minTRGG, s0, 'OutputForm', 'cell');
    
    disp(['...[' num2str(idSS) '] vector representation for each vertex']);
    
    % ---------------
    % ===For TREE===
    % vector representation for each vertex 1 --> nGG

    disp('......vector representation for each vertex');
    % length(wwGG): #edges in graph GG (can be reduced into #edges in tree)

    vecGG_VV = zeros(nGG, length(minTR_EdgeWW));
    
    for ii = 1:nGG % each vertex in graph
        vecGG_VV(ii, trEP{ii}) = 1;
    end

    % ALREADY tree!!!

    disp('......vector representation for each distribution');
    % ===For Data===
    % N: #samples (input data)
    % Input: WW, 

    % V2: --> spare version
    XX_SI = zeros(N, length(minTR_EdgeWW));
        
    % sum mass
    XX_mass = zeros(N, 1); % column vector
    
    for ii = 1:N % each distribution
        
        % == For unbalanced version (without normalization)
        tmpWW = WW{ii}; % without normalization for weight!!!
        XX_mass(ii) = sum(WW{ii});
        
        tmpXX = XX_ID{ii}; % id of vertices (multiset due to random graph)
        
        tmpXX_GG_TR = vecGG_VV(tmpXX, :);
        
        tmpWW_GG_TR = repmat(tmpWW, 1, length(minTR_EdgeWW));
        
        tmpWWXX = tmpXX_GG_TR .* tmpWW_GG_TR;
        XX_SI(ii, :) = sum(tmpWWXX, 1);
    end
    runTime_Prep(idSS) = toc;
   
    tic
    % compute the Lp distance matrix
    DD_SS_II = zeros(N, N); % unbalanced Sobolev
    varDD_SS_II = zeros(N, N); % variant of unbalanced Sobolev (minus summarization term)
    
    for ii = 1:N
        
        % % ii --> (ii+1):N
        % ii --> ii:N
        if mod(ii, 20) == 0
            disp(['...' num2str(ii)]);
        end
    
        tmpII_vec = XX_SI(ii, :);
        
        tmpJJ_mat = XX_SI(ii:N, :);
        tmpII_mat = repmat(tmpII_vec, N-ii+1, 1);

        tmpAbsDD_mat = abs(tmpII_mat - tmpJJ_mat);
        
        if pp > 1
            tmpPP_AbsDD_mat = tmpAbsDD_mat.^pp;
        else
            tmpPP_AbsDD_mat = tmpAbsDD_mat;
        end
        
        wwGG_TR_mat = repmat(minTR_EdgeWW', N-ii+1, 1); 
        
        % --
        tmpWWPP_AbsDD_mat = wwGG_TR_mat .* tmpPP_AbsDD_mat;
        
        tmpPP_DD_vec = sum(tmpWWPP_AbsDD_mat, 2); % sum over rows --> column
        
        % p = 1
        tmpDD_vec = tmpPP_DD_vec;

        % for the summarization over tree
        % tmpDD_vec: column vector!!!
        
        %------------
        wwSub = opt.a + (opt.b*opt.lambda)*0.5 - opt.alpha;
        tmpII_Sub = repmat(XX_mass(ii), N-ii+1, 1);
        tmpJJ_Sub = XX_mass(ii:N);
        
        tmpSub_vec = wwSub*abs(tmpII_Sub - tmpJJ_Sub); % column vector
        
        %------------
        wwSum = opt.b*opt.lambda*0.5;
        tmpSum_vec = wwSum*(tmpII_Sub + tmpJJ_Sub); % column vector

        %------------
        tmpDD_vec_all = opt.b*tmpDD_vec + tmpSub_vec;
        vartmpDD_vec_all = tmpDD_vec_all - tmpSum_vec;
         
        DD_SS_II(ii, ii:N) = tmpDD_vec_all';
        DD_SS_II(ii:N, ii) = tmpDD_vec_all;
        
        varDD_SS_II(ii, ii:N) = vartmpDD_vec_all';
        varDD_SS_II(ii:N, ii) = vartmpDD_vec_all;
        
    end
    runTime_Dist(idSS) = toc;
    
    % save distance matrix
    DD_SS{idSS} = DD_SS_II;
    varDD_SS{idSS} = varDD_SS_II;
end

runTime_TreeSampling_Avg = sum(runTime_TreeSampling) / nSS;
runTime_Prep_Avg = sum(runTime_Prep) / nSS;
runTime_Dist_Avg = sum(runTime_Dist) / nSS;

runTime_Dist_ALL = runTime_Prep + runTime_Dist + runTime_TreeSampling;
runTime_Dist_ALL_Avg = sum(runTime_Dist_ALL) / nSS;

% Average
tmpNN = [1, 5, 10, 20];
tmpDDSS_Cell = cell(length(tmpNN), 1);
vartmpDDSS_Cell = cell(length(tmpNN), 1);

for iiRR = 1:length(tmpNN)
    
    tmpDDSS = zeros(N, N);
    vartmpDDSS = zeros(N, N);

    for ii = 1:tmpNN(iiRR)
        tmpDDSS = tmpDDSS + DD_SS{ii};
        vartmpDDSS = vartmpDDSS + varDD_SS{ii};
    end
    
    tmpDDSS = tmpDDSS / tmpNN(iiRR);
    vartmpDDSS = vartmpDDSS / tmpNN(iiRR);
    
    tmpDDSS_Cell{iiRR} = tmpDDSS;
    vartmpDDSS_Cell{iiRR} = vartmpDDSS;
end

% ---> dAlpha
DD_SS1 = tmpDDSS_Cell{1};
DD_SS5 = tmpDDSS_Cell{2};
DD_SS10 = tmpDDSS_Cell{3};
DD_SS20 = tmpDDSS_Cell{4};

%-------
% ---> EPT
varDD_SS1 = vartmpDDSS_Cell{1};
varDD_SS5 = vartmpDDSS_Cell{2};
varDD_SS10 = vartmpDDSS_Cell{3};
varDD_SS20 = vartmpDDSS_Cell{4};

outName = [dsName '_RandTree_dAlpha_EPT_' num2str(maxKC) '_' typeGG '_4S_P' num2str(pp) '.mat'];

save(outName, 'DD_SS1', 'DD_SS5', 'DD_SS10', 'DD_SS20', ...
     'varDD_SS1', 'varDD_SS5', 'varDD_SS10', 'varDD_SS20', ...
     'runTime_TreeSampling', 'runTime_Dist', 'runTime_Prep', 'runTime_Dist_ALL', ...
     'runTime_TreeSampling_Avg', 'runTime_Dist_Avg', 'runTime_Prep_Avg', 'runTime_Dist_ALL_Avg', ...
     'randSArray', 'pp', 'nSS', ...
     'YY');

disp('======================================');

end
end

disp('FINISH ALL !!!');

