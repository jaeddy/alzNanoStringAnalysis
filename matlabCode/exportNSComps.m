% Export DEG results for IPA

[comps{1:numel(degStruct)}] = deal(degStruct.comp);
comps = comps';

threshold = 0.05;

ipaMat = repmat({''},numel(names)+1,numel(degStruct));
for i = 1:numel(degStruct)
    entry = [degStruct(i).comp; degStruct(i).results(:,3)];
    if threshold > 0
        qvals = cellfun(@str2num, degStruct(i).results(:,4));
        qvals(isnan(qvals)) = 1;
        entry(find(qvals > threshold)+1) = {'0'};
    end
    ipaMat(:,i) = entry;
end
ipaMat = [[{'Names'};names], ipaMat];


ipaMatSmall = ipaMat(:,[1,2,14,2,3]);
ipaMatSmall(1,[2,4]) = {'Tg_IL10 vs. Tg_CTRL'};
ipaMatSmall(1,[3]) = {'WT_IL10 vs. WT_CTRL'};
ipaMatSmall(1,[5]) = {'Tg_OLD vs. Tg_CTRL'};

%%
ipaMatSmall(~ismember(ipaMatSmall(:,3),'0'),2) = {'0'};
ipaMatSmall(1,2) = {'(Tg_IL10vs.Tg_CTRL)-(WT_IL10vs.WT_CTRL)'};

ipaMatSmall(~ismember(ipaMatSmall(:,5),'0'),4) = {'0'};
ipaMatSmall(1,4) = {'(Tg_IL10vs.Tg_CTRL)-(Tg_OLDvs.Tg_CTRL)'};

ipaMatSmall(:,[3,5]) = [];
for i = 2:3
    ipaMatSmall(2:end,i) = cellstr(num2str(...
        -1.*cellfun(@str2num,ipaMatSmall(2:end,i))));
end

textArrayWrite('subComps_qvalThres.txt',ipaMatSmall,'\t')

% if threshold > 0
%     textArrayWrite('allComps_qvalThres.txt',ipaMat,'\t')
% else textArrayWrite('allComps_noThres.txt',ipaMat,'\t')
% end

%% compile comp results for Excel

headers = {'Accession','IPA Name','Log ratio','P-value','Q-value','FDR'};
for i = 1:numel(degStruct)
    degTable = degStruct(i).results;
    degAccs = degTable(:,1);
    degNames = name2accession(...
        directIntersect(name2accession(:,2),degAccs),1);
    degTable = [degTable(:,1),degNames,degTable(:,2:end)];
    degTable = [headers; degTable];
    xlsStruct(i).comp = degStruct(i).comp;
    xlsStruct(i).table = degTable;
end

%% subtractive comps with DIRAC

% subtractive comp 1, VE
TgIL10vsTgCTRL = diracStruct(1).results_tables.Variably_Expressed_Networks;
WTIL10vsWTCTRL = diracStruct(2).results_tables.Variably_Expressed_Networks;
TgIL10vsTgCTRL(90:end,:) = [];
WTIL10vsWTCTRL(98:end,:) = [];
TgIL10vsTgCTRL(...
    find(ismember(TgIL10vsTgCTRL(:,1),WTIL10vsWTCTRL(:,1)))+1,:) = [];

% subtractive comp 1, DR
TgIL10vsTgCTRL = ...
    diracStruct(1).results_tables.Differentially_Regulated_Networks;
WTIL10vsWTCTRL = ...
    diracStruct(2).results_tables.Differentially_Regulated_Networks;
TgIL10vsTgCTRL(28:end,:) = [];
WTIL10vsWTCTRL(24:end,:) = [];
TgIL10vsTgCTRL(...
    find(ismember(TgIL10vsTgCTRL(:,1),WTIL10vsWTCTRL(:,1)))+1,:) = [];

% subtractive comp 2, VE
TgIL10vsTgCTRL = diracStruct(1).results_tables.Variably_Expressed_Networks;
TgOLDvsTgCTRL = diracStruct(4).results_tables.Variably_Expressed_Networks;
TgIL10vsTgCTRL(90:end,:) = [];
TgOLDvsTgCTRL(71:end,:) = [];
TgIL10vsTgCTRL(...
    find(ismember(TgIL10vsTgCTRL(:,1),TgOLDvsTgCTRL(:,1)))+1,:) = [];

% subtractive comp 1, DR
TgIL10vsTgCTRL = ...
    diracStruct(1).results_tables.Differentially_Regulated_Networks;
TgOLDvsTgCTRL = ...
    diracStruct(4).results_tables.Differentially_Regulated_Networks;
TgIL10vsTgCTRL(28:end,:) = [];
TgOLDvsTgCTRL(2:end,:) = [];
TgIL10vsTgCTRL(...
    find(ismember(TgIL10vsTgCTRL(:,1),TgOLDvsTgCTRL(:,1)))+1,:) = [];