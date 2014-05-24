function degStruct = degAllComps(X,Y,names)

% This function is designed to perform differentially expressed gene
% calculations on all pairs of sample groups. A two-sample t-test is used
% in the mattest function to obtain p-values for all genes, and the mafdr
% function is used to estimate multiple test correction FDR and q-values.

sampleSubsets = unique(Y);
sampleComps = nchoosek(1:numel(sampleSubsets),2);

for i = 1:size(sampleComps,1)
    group1 = strcmp(sampleSubsets{sampleComps(i,1)},Y);
    group2 = strcmp(sampleSubsets{sampleComps(i,2)},Y);

    pvals = mattest(X(:,group1),X(:,group2));
    [~,sortIdx] = sort(pvals,'ascend');
%     sortIdx = 1:numel(pvals);

    foldChange = mean(X(:,group1),2) - mean(X(:,group2),2);
        
    [fdr,qvals] = mafdr(pvals);

    degTable = [names, cellstr(num2str(pvals)),...
        cellstr(num2str(foldChange)),...
        cellstr(num2str(qvals)), cellstr(num2str(fdr))];
    degTable = degTable(sortIdx,:);
    
    degStruct(i).comp = [sampleSubsets{sampleComps(i,1)},' vs. ',...
        sampleSubsets{sampleComps(i,2)}];
    degStruct(i).results = degTable;
end