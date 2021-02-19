% 02/08/2021 Makoto. Finished. Shaun's Excel data are screwed up. Fixed.
% 02/05/2021 Makoto. Created.

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


wHandle = waitbar(0, 'Starting...');
icCounter = 0;
for setIdx = 1:length(subSetIdx)

    waitbar(setIdx/length(subSetIdx), wHandle)
    
    % 1. Load data.
    loadName = allSets(subSetIdx(setIdx)).name;
    EEG = pop_loadset('filename', loadName, 'filepath', '/data/mobi/Hiroki/p4000_import_preprocess_n600', 'loadmode', 'info');

    % 2. Obtain datapoint rejection rate by ASR.
    subjectwiseCleaningStats(setIdx,1) = sum(EEG.etc.clean_sample_mask)/length(EEG.etc.clean_sample_mask);

    % 3. Obtain the number of ICs left for three levels of cleaning.
    subjectwiseCleaningStats(setIdx,2) = size(EEG.icaweights,1);
    subjectwiseCleaningStats(setIdx,3) = size(EEG.Hiroki.icRej1.icaweights,1);
    subjectwiseCleaningStats(setIdx,4) = size(EEG.Hiroki.icRej2.icaweights,1);

    % 4. Obtain IC-wise properties.
    load(EEG.dipfit.hdmfile); % This returns 'vol'.
    dipoleXyz = zeros(length(EEG.dipfit.model),3);
    for icIdx = 1:length(EEG.dipfit.model)
        icCounter = icCounter + 1;
        dipoleXyz = EEG.dipfit.model(icIdx).posxyz(1,:);
        depth     = ft_sourcedepth(dipoleXyz, vol);
        inputData = [setIdx, icIdx, dipoleXyz, size(EEG.dipfit.model(icIdx).posxyz,1), EEG.dipfit.model(icIdx).rv, depth, EEG.etc.ic_classification.ICLabel.classifications(icIdx,:)];
        icwiseData(icCounter,:) = inputData;
    end

    % 5. Obtain PSD.
    psd_orig(setIdx,:,:)   = EEG.Hiroki.precleanedEEG.spectra;
    psd_asr(setIdx,:,:)    = EEG.Hiroki.asrCleaned.spectra;
    psd_clean1(setIdx,:,:) = EEG.Hiroki.icRej1.spectra;
    psd_clean2(setIdx,:,:) = EEG.Hiroki.icRej2.spectra;
    
    % 6. Obtain MI and resultant vector length.
    pacTensor(setIdx,:,1) = EEG.Hiroki.precleanedEEG.pac.mi;
    pacTensor(setIdx,:,2) = EEG.Hiroki.icRej1.pac.mi;
    pacTensor(setIdx,:,3) = EEG.Hiroki.icRej2.pac.mi;
    pacTensor(setIdx,:,4) = EEG.Hiroki.precleanedEEG.pac.vectLength;
    pacTensor(setIdx,:,5) = EEG.Hiroki.icRej1.pac.vectLength;
    pacTensor(setIdx,:,6) = EEG.Hiroki.icRej2.pac.vectLength;
    
    % 7. Obtain group index.
    currentSubjIdx = subSetIdx(setIdx);
    if any(intersect(preWake1MatchedSubjIdx,currentSubjIdx))
        pacGroupIdx(setIdx) = 1;
    elseif any(intersect(preSleep1MatchedSubjIdx,currentSubjIdx))
        pacGroupIdx(setIdx) = 2;
    elseif any(intersect(preWake2MatchedSubjIdx,currentSubjIdx))
        pacGroupIdx(setIdx) = 3;
    elseif any(intersect(preSleep2MatchedSubjIdx,currentSubjIdx))
        pacGroupIdx(setIdx) = 4;
    end
end
save /data/mobi/Hiroki/p4200_select50subjects/summaryData subjectwiseCleaningStats icwiseData psd_orig psd_asr psd_clean1 psd_clean2 pacTensor pacTensor pacGroupIdx