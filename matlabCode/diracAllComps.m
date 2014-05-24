function diracStruct = diracAllComps(X,Y,names)

% This function is designed to perform differentially expressed gene
% calculations on all pairs of sample groups. A two-sample t-test is used
% in the mattest function to obtain p-values for all genes, and the mafdr
% function is used to estimate multiple test correction FDR and q-values.

geneset_defs_file = ...
    ['/Users/jameslocal/Google Drive/MyCode/MatlabCode/',...
    'Expression Analysis/AD Expression/geneset_definitions_2014.mat'];
% geneset_defs_opt = 'canonpath2014_gs_defs';
geneset_defs_opt = 'ipaCanonPaths';

gs_min = 5;
top_gs_M = [];
num_permutations = 1000;

sampleSubsets = unique(Y);
sampleComps = nchoosek(1:numel(sampleSubsets),2);

sampleComps = sampleComps([1,13,3,2,14,6],:);
sampleComps([1,2,4,5],:) = fliplr(sampleComps([1,2,4,5],:));

for i = 1:size(sampleComps,1)
    group1 = strcmp(sampleSubsets{sampleComps(i,1)},Y);
    group2 = strcmp(sampleSubsets{sampleComps(i,2)},Y);
    
    X_i = [X(:,group1),X(:,group2)];
    Y_i = [Y(group1),Y(:,group2)];
    
    groups = strcmp(sampleSubsets{sampleComps(i,1)},Y_i);

    [results_tables,results_raw] = ...
        dirac(X_i,names,groups,geneset_defs_file,geneset_defs_opt,...
        gs_min,top_gs_M,num_permutations);
    
    diracStruct(i).comp = [sampleSubsets{sampleComps(i,1)},'vs.',...
        sampleSubsets{sampleComps(i,2)}];
    diracStruct(i).results_tables = results_tables;
    diracStruct(i).results_raw = results_raw;
end

save('diracResults.mat','diracStruct')