% 02/09/2021 Makoto. created.

clear

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract 50 subjects for each of 2x2 factorial design. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NUM,TXT,RAW] = xlsread('/data/mobi/Hiroki/code/QHYPS_020221_1647_exportMM.xlsx');

preWake1Idx  = find(contains(TXT(1,:),'CodePreWake1'));
preSleep1Idx = find(contains(TXT(1,:),'CodePreSleep1'));
preWake2Idx  = find(contains(TXT(1,:),'CodePreWake2'));
preSleep2Idx = find(contains(TXT(1,:),'CodePreSleep2'));

preWake1ExcelSubjNames  = sort(TXT(2:101,preWake1Idx));
preSleep1ExcelSubjNames = sort(TXT(2:101,preSleep1Idx));
preWake2ExcelSubjNames  = sort(TXT(2:101,preWake2Idx));
preSleep2ExcelSubjNames = sort(TXT(2:101,preSleep2Idx));

allSets = dir('/data/mobi/Hiroki/p4000_import_preprocess_n600/*.set');
allSetNames = {allSets.name}';
allSetNames = cellfun(@(x) x(1:end-4), allSetNames, 'uniformoutput', false);

preWake1MatchedSubjIdx  = find(contains(allSetNames, preWake1ExcelSubjNames));
preSleep1MatchedSubjIdx = find(contains(allSetNames, preSleep1ExcelSubjNames));
preWake2MatchedSubjIdx  = find(contains(allSetNames, preWake2ExcelSubjNames));
preSleep2MatchedSubjIdx = find(contains(allSetNames, preSleep2ExcelSubjNames));

preWake1MatchedSubjNames  = allSetNames(preWake1MatchedSubjIdx);
preSleep1MatchedSubjNames = allSetNames(preSleep1MatchedSubjIdx);
preWake2MatchedSubjNames  = allSetNames(preWake2MatchedSubjIdx);
preSleep2MatchedSubjNames = allSetNames(preSleep2MatchedSubjIdx);

noExistingWake1Names  = setdiff(preWake1ExcelSubjNames,  preWake1MatchedSubjNames);
noExistingSleep1Names = setdiff(preSleep1ExcelSubjNames, preSleep1MatchedSubjNames);
noExistingWake2Names  = setdiff(preWake2ExcelSubjNames,  preWake2MatchedSubjNames);
noExistingSleep2Names = setdiff(preSleep2ExcelSubjNames, preSleep2MatchedSubjNames);

% Names in {} are from the Excel file.
% % Wake 1
% {'DQB93a'} -> DQB93 exists.
% {'EFV94' } -> efv94 exists.
% {'QGC43a'} -> QGC43 exists.
% {'TWG47' } -> twg47 exists.
% {'XBI94a'} -> XBI94 exists.
% 
% % Sleep 1
% {'JDG25a'} -> JDG25 exists.
% {'STV75a'} -> STV75 exists.
% 
% % Wake 2
% {'EGD95a'} -> EGD95 exists.
% {'USP38a'} -> USP38 exists.
% 
% % Sleep 2
% {'DEV35a'} -> DEV35 exists.
% {'TDS79a'} -> TDS79 exists.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Rename the Excel-derived name list to match the set files. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
renamingIdx = find(contains(preWake1ExcelSubjNames, 'DQB93a'));  preWake1ExcelSubjNames{renamingIdx} = 'DQB93';
renamingIdx = find(contains(preWake1ExcelSubjNames, 'EFV94'));   preWake1ExcelSubjNames{renamingIdx} = 'efv94';
renamingIdx = find(contains(preWake1ExcelSubjNames, 'QGC43a'));  preWake1ExcelSubjNames{renamingIdx} = 'QGC43';
renamingIdx = find(contains(preWake1ExcelSubjNames, 'TWG47'));   preWake1ExcelSubjNames{renamingIdx} = 'twg47';
renamingIdx = find(contains(preWake1ExcelSubjNames, 'XBI94a'));  preWake1ExcelSubjNames{renamingIdx} = 'XBI94';

renamingIdx = find(contains(preSleep1ExcelSubjNames, 'JDG25a')); preSleep1ExcelSubjNames{renamingIdx} = 'JDG25';
renamingIdx = find(contains(preSleep1ExcelSubjNames, 'STV75a')); preSleep1ExcelSubjNames{renamingIdx} = 'STV75';

renamingIdx = find(contains(preWake2ExcelSubjNames, 'EGD95a'));  preWake2ExcelSubjNames{renamingIdx} = 'EGD95';
renamingIdx = find(contains(preWake2ExcelSubjNames, 'USP38a'));  preWake2ExcelSubjNames{renamingIdx} = 'USP38';

renamingIdx = find(contains(preSleep2ExcelSubjNames, 'DEV35a')); preSleep2ExcelSubjNames{renamingIdx} = 'DEV35';
renamingIdx = find(contains(preSleep2ExcelSubjNames, 'TDS79a')); preSleep2ExcelSubjNames{renamingIdx} = 'TDS79';

% Re-define the dataset indices.
preWake1MatchedSubjIdx  = find(contains(allSetNames, preWake1ExcelSubjNames));
preSleep1MatchedSubjIdx = find(contains(allSetNames, preSleep1ExcelSubjNames));
preWake2MatchedSubjIdx  = find(contains(allSetNames, preWake2ExcelSubjNames));
preSleep2MatchedSubjIdx = find(contains(allSetNames, preSleep2ExcelSubjNames));

subSetIdx = sort([preWake1MatchedSubjIdx; preSleep1MatchedSubjIdx; preWake2MatchedSubjIdx; preSleep2MatchedSubjIdx]);

subjectwiseCleaningStats = zeros(length(subSetIdx),4); % dataRejectionRate, numICs_cleaning0-1-2.
icwiseData = zeros(1,15); % subjIdx, icIdx, x, y, z, singleOrDouble, rv, dipoleDepth, classLabelProbability
psd_orig   = zeros(length(subSetIdx), 21, 401);
psd_asr    = zeros(length(subSetIdx), 21, 401);
psd_clean1 = zeros(length(subSetIdx), 21, 401);
psd_clean2 = zeros(length(subSetIdx), 21, 401);

pacTensor   = zeros(400, 21, 6);
pacGroupIdx = zeros(400, 1);




allSets = dir('/data/mobi/Hiroki/p4000_import_preprocess_n600/*.set');

varRatioTensor = zeros(21,6,length(subSetIdx));
for setIdx = 1:length(subSetIdx)
    
    % 1. Load data.
    loadName = allSets(subSetIdx(setIdx)).name;
    EEG = pop_loadset('filename', loadName, 'filepath', '/data/mobi/Hiroki/p4000_import_preprocess_n600');
    
    % 2. Calculate variance difference by ASR alone.
    asrMask = EEG.etc.clean_sample_mask;
    rawVar    = var(EEG.Hiroki.precleanedEEG.data(:,asrMask), 0, 2);
    asrVar    = var(EEG.data, 0, 2);
    icrej1Var = var(EEG.Hiroki.icRej1.data, 0, 2);
    icrej2Var = var(EEG.Hiroki.icRej2.data, 0, 2);
    
    varRatioTensor(:,1,setIdx) = asrVar./rawVar;
    varRatioTensor(:,2,setIdx) = icrej1Var./rawVar;
    varRatioTensor(:,3,setIdx) = icrej2Var./rawVar;
    varRatioTensor(:,4,setIdx) = icrej1Var./asrVar;
    varRatioTensor(:,5,setIdx) = icrej2Var./asrVar;
    varRatioTensor(:,6,setIdx) = icrej2Var./icrej1Var;
end
save /data/mobi/Hiroki/p4220_calculateVarianceChange/varRatioTensor varRatioTensor

%%

asr_Raw = (1-mean(varRatioTensor(:,1,:),3))*-1;
ic1_Raw = (1-mean(varRatioTensor(:,2,:),3))*-1;
ic2_Raw = (1-mean(varRatioTensor(:,3,:),3))*-1;
ic1_asr = (1-mean(varRatioTensor(:,4,:),3))*-1;
ic2_asr = (1-mean(varRatioTensor(:,5,:),3))*-1;
ic2_ic1 = (1-mean(varRatioTensor(:,6,:),3))*-1;

asr_Raw_std = std(varRatioTensor(:,1,:),0,3);
ic1_Raw_std = std(varRatioTensor(:,2,:),0,3);
ic2_Raw_std = std(varRatioTensor(:,3,:),0,3);
ic1_asr_std = std(varRatioTensor(:,4,:),0,3);
ic2_asr_std = std(varRatioTensor(:,5,:),0,3);
ic2_ic1_std = std(varRatioTensor(:,6,:),0,3);

    % figure
    % subplot(3,3,1)
    % topoplot(asr_Raw, EEG.chanlocs, 'maplimits', [-0.8 0.8])
    % subplot(3,3,2)
    % topoplot(ic1_Raw, EEG.chanlocs, 'maplimits', [-0.8 0.8])
    % subplot(3,3,3)
    % topoplot(ic2_Raw, EEG.chanlocs, 'maplimits', [-0.8 0.8])
    % 
    % subplot(3,3,5)
    % topoplot(ic1_asr, EEG.chanlocs, 'maplimits', [-0.8 0.8])
    % subplot(3,3,6)
    % topoplot(ic2_asr, EEG.chanlocs, 'maplimits', [-0.8 0.8])
    % 
    % subplot(3,3,9)
    % topoplot(ic2_ic1, EEG.chanlocs, 'maplimits', [-0.8 0.8])
    % 
    % currentPosition = get(gca, 'position')
    % colorbar
    % set(gca, 'position', currentPosition)

% Z-score version
figure
subplot(3,3,1)
topoplot(zscore(asr_Raw), EEG.chanlocs, 'maplimits', [-3 3])
title(sprintf('var. reduc.\n%.0f-%.0f%% (SD %.0f)', abs(max(asr_Raw)*100), abs(min(asr_Raw)*100), mean(asr_Raw_std)*100))

subplot(3,3,2)
topoplot(zscore(ic1_Raw), EEG.chanlocs, 'maplimits', [-3 3])
title(sprintf('var. reduc.\n%.0f-%.0f%% (SD %.0f)', abs(max(ic1_Raw)*100), abs(min(ic1_Raw)*100), mean(ic1_Raw_std)*100))

subplot(3,3,3)
topoplot(zscore(ic2_Raw), EEG.chanlocs, 'maplimits', [-3 3])
title(sprintf('var. reduc.\n%.0f-%.0f%% (SD %.0f)', abs(max(ic2_Raw)*100), abs(min(ic2_Raw)*100), mean(ic2_Raw_std)*100))

subplot(3,3,5)
topoplot(zscore(ic1_asr), EEG.chanlocs, 'maplimits', [-3 3])
title(sprintf('var. reduc.\n%.0f-%.0f%% (SD %.0f)', abs(max(ic1_asr)*100), abs(min(ic1_asr)*100), mean(ic1_asr_std)*100))

subplot(3,3,6)
topoplot(zscore(ic2_asr), EEG.chanlocs, 'maplimits', [-3 3])
title(sprintf('var. reduc.\n%.0f-%.0f%% (SD %.0f)', abs(max(ic2_asr)*100), abs(min(ic2_asr)*100), mean(ic2_asr_std)*100))

subplot(3,3,9)
topoplot(zscore(ic2_ic1), EEG.chanlocs, 'maplimits', [-3 3])
title(sprintf('var. reduc.\n%.0f-%.0f%% (SD %.0f)', abs(max(ic2_ic1)*100), abs(min(ic2_ic1)*100), mean(ic2_ic1_std)*100))

currentPosition = get(gca, 'position')
colorbar
set(gca, 'position', currentPosition)

% Print.
set(gcf, 'position', [1           2        1858         929])
print('/data/mobi/Hiroki/p4220_calculateVarianceChange/varChanges', '-dsvg')