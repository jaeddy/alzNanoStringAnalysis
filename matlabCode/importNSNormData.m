% Import normalized AD NanoString data

%% Read text file output from NanoStringNorm
fileDir = '/Users/jeddy/Google Drive/MyCode/RCode/AD_nanoString/results/';
fileName = 'AD_NanoStringNorm_data_hkLog.txt';

numCol = 31;
delim = '\t';

AD_NS_norm = textArrayRead([fileDir,fileName],numCol,delim);

%% Clean up formatting in main table
AD_NS_norm = regexprep(AD_NS_norm,'"','');
AD_NS_norm = strtrim(AD_NS_norm);
AD_NS_norm(2:end,1:end-1) = AD_NS_norm(2:end,2:end);
AD_NS_norm(:,end) = [];

%% Pull out and format variables for data, row names, class labels
X = cellfun(@str2num, AD_NS_norm(2:end,4:end));
names = AD_NS_norm(2:end,2);

% Define class labels
Y = repmat({''},1,27);
Y([1:3,7:9]) = {'Tg_IL10'};
Y([4:6,10:12]) = {'Tg_CTRL'};
Y([13,17:19]) = {'WT_IL10'};
Y([14:16,20:21]) = {'WT_CTRL'};
Y([22:24]) = {'Tg_OLD'};
Y([25:27]) = {'WT_OLD'};

%% Remove control probes
negControls = ~cellfun('isempty',strfind(names,'NEG'));
posControls = ~cellfun('isempty',strfind(names,'POS'));
allControls = negControls | posControls;

X(allControls,:) = [];
names(allControls) = [];

%% Save data
save AD_NSNorm_data X Y names
