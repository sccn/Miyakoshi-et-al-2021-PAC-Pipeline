% 01/25/2021 Makoto. Used.
% 01/19/2021 Makoto. Used.
% 11/13/2020 Makoto. Created.

% addpath('/data/projects/makoto/Tools/mdiToolbox4')

allSets = dir('/data/mobi/Hiroki/p4000_import_preprocess_n600/*.set');
subjNameList = cell(1,length(allSets));
pacRaw    = nan(21,length(allSets));
pacClean1 = nan(21,length(allSets));
pacClean2 = nan(21,length(allSets));
psdRaw    = nan(21,401,length(allSets));
psdClean0 = nan(21,401,length(allSets));
psdClean1 = nan(21,401,length(allSets));
psdClean2 = nan(21,401,length(allSets));
numIcs1 = nan(length(allSets),1);
numIcs2 = nan(length(allSets),1);
dataLength = nan(length(allSets),1);
for setIdx = 1:length(allSets)
    
    loadName = allSets(setIdx).name;
    dataName = loadName(1:end-4);
    EEG      = pop_loadset('filename', loadName, 'filepath', '/data/mobi/Hiroki/p4000_import_preprocess_n600', 'loadmode', 'info');
    
    % Store final output.
    subjNameList{1,setIdx} = dataName;
    pacRaw(:,setIdx)    = EEG.Hiroki.precleanedEEG.pac.mi;
    pacClean1(:,setIdx) = EEG.Hiroki.icRej1.pac.mi;
    pacClean2(:,setIdx) = EEG.Hiroki.icRej2.pac.mi;
    
    psdRaw(:,:,setIdx)    = EEG.Hiroki.precleanedEEG.spectra;
    psdClean0(:,:,setIdx) = EEG.Hiroki.asrCleaned.spectra;
    psdClean1(:,:,setIdx) = EEG.Hiroki.icRej1.spectra;
    psdClean2(:,:,setIdx) = EEG.Hiroki.icRej2.spectra;
    
    numIcs1(setIdx) = size(EEG.Hiroki.icRej1.icaweights,1);
    numIcs2(setIdx) = size(EEG.Hiroki.icRej2.icaweights,1);
    dataLength(setIdx) = EEG.xmax;
end

save /data/mobi/Hiroki/p4100_groupLevelStructure/fullChData.mat subjNameList pacRaw pacClean1 pacClean2 psdRaw psdClean0 psdClean1 psdClean2 numIcs1 numIcs2 dataLength
       
%% Save the imputed data.
chanLabels = {EEG.chanlocs.labels}';
chanLabels{end+1,1} = 'Global MI';

rawDataTable = array2table([pacRaw; sum(pacRaw,1)],...
                           'VariableNames', subjNameList,...
                           'RowNames', chanLabels);
clean1DataTable = array2table([pacClean1; sum(pacClean1,1)],...
                           'VariableNames', subjNameList,...
                           'RowNames', chanLabels);
clean2DataTable = array2table([pacClean2; sum(pacClean2,1)],...
                           'VariableNames', subjNameList,...
                           'RowNames', chanLabels);

writetable(rawDataTable,   '/data/mobi/Hiroki/p4100_groupLevelStructure/rawData.xlsx',...
                           'FileType', 'spreadsheet',...
                           'WriteVariableNames', true,...
                           'WriteRowNames', true);
writetable(clean1DataTable, '/data/mobi/Hiroki/p4100_groupLevelStructure/clean1Data.xlsx',...
                           'FileType', 'spreadsheet',...
                           'WriteVariableNames', true,...
                           'WriteRowNames', true);
writetable(clean2DataTable, '/data/mobi/Hiroki/p4100_groupLevelStructure/clean2Data.xlsx',...
                           'FileType', 'spreadsheet',...
                           'WriteVariableNames', true,...
                           'WriteRowNames', true);

                       
                       
%%
allSets = dir('/data/mobi/Hiroki/p4000_import_preprocess_n600/*.set');
subjectwiseCleaningStats = zeros(length(allSets),4); % dataRejectionRate, numICs_cleaning0-1-2.
icwiseData = zeros(1,15); % subjIdx, icIdx, x, y, z, singleOrDouble, rv, dipoleDepth, classLabelProbability
psd_orig   = zeros(length(allSets), 21, 401);
psd_asr    = zeros(length(allSets), 21, 401);
psd_clean1 = zeros(length(allSets), 21, 401);
psd_clean2 = zeros(length(allSets), 21, 401);

icCounter = 0;
for setIdx = 1:length(allSets)
    
    % 1. Load data.
    loadName = allSets(setIdx).name;
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
end

save /data/mobi/Hiroki/p4100_groupLevelStructure icwiseData psd_orig psd_asr psd_clean1 psd_clean2



%% Generate figures

load /data/mobi/Hiroki/p4100_groupLevelStructure/groupLevelStructure.mat
addpath /data/projects/makoto/Tools/RainCloudPlots-master/tutorial_matlab
addpath /data/mobi/Daisuke/code


% 1. Histogram of ASR datapoint rejections.
figure
histHandle = histogram([1-subjectwiseCleaningStats(:,1)]*100, 0:5:50, 'FaceColor', [0.66 0.76 1]);
title(sprintf('ASR rejection rate (M=%.1f, SD=%.1f, range %.0f-%.1f)',...
              mean(1-subjectwiseCleaningStats(:,1))*100, ...
              std(1-subjectwiseCleaningStats(:,1))*100, ...
              min(1-subjectwiseCleaningStats(:,1))*100, ...
              max(1-subjectwiseCleaningStats(:,1))*100))
xlabel('Datapoint rejection rate (%)')
ylabel('Number of datasets')
ylim([0 395])
axis square
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 24)
set(gcf, 'position', [1           2        1858         929])
%print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure1.jpg', '-r300', '-djpeg95')
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure2', '-dsvg')



% 2. Number of ICs rejected.
figure
subplot(1,2,1)
histHandle = histogram(21-subjectwiseCleaningStats(:,3), 0:1:21, 'FaceColor', [0.66 0.76 1]);
title(sprintf('Cleaning Level 1:\nNumber of ICs rejected\nM=%.1f (SD=%.1f, range %.0f-%.1f)',...
              mean(21-subjectwiseCleaningStats(:,3)), ...
              std(21-subjectwiseCleaningStats(:,3)), ...
              min(21-subjectwiseCleaningStats(:,3)), ...
              max(21-subjectwiseCleaningStats(:,3))))
xlabel('Number of ICs rejected')
ylabel('Number of datasets')
axis square
xlim([1 21])

subplot(1,2,2)
histHandle = histogram(21-subjectwiseCleaningStats(:,4), 0:1:21, 'FaceColor', [0.66 0.76 1]);
title(sprintf('Cleaning Level 2:\nNumber of ICs rejected\nM=%.1f (SD=%.1f, range %.0f-%.1f))',...
              mean(21-subjectwiseCleaningStats(:,4)), ...
              std(21-subjectwiseCleaningStats(:,4)), ...
              min(21-subjectwiseCleaningStats(:,4)), ...
              max(21-subjectwiseCleaningStats(:,4))))
xlabel('Number of ICs')
ylabel('Number of datasets')
axis square
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 24)
set(gcf, 'position', [1           2        1858         929])
xlim([1 21])
% print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure2.jpg', '-r300', '-djpeg95')
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure3', '-dsvg')



% 3. Class label rate. subjIdx, icIdx, x, y, z, singleOrDouble, rv, dipoleDepth, classLabelProbability
[~,topLabelIdx] = max(icwiseData(:,9:end), [], 2);
classLabelCounts =  [sum(topLabelIdx==1) ...
                     sum(topLabelIdx==2) ...
                     sum(topLabelIdx==3) ...
                     sum(topLabelIdx==4) ...
                     sum(topLabelIdx==5) ...
                     sum(topLabelIdx==6) ...
                     sum(topLabelIdx==7)];
                 
classLabelCountsSorted = sort(classLabelCounts, 'descend');
classLabelCountsSorted([2:7]) = classLabelCountsSorted([3 4 5 6 7 2]);
                 
figure
explode = [1 0 0 0 0 0 1];
pieHandle = pie(classLabelCountsSorted, explode);

% Obtain fancy colors from cbrewer.
addpath('/data/projects/makoto/Tools/cbrewer/cbrewer')
colormap = cbrewer('qual', 'Set3', 7);

colormap([2 7], :) = colormap([7 2], :);
colormap([2 6], :) = colormap([6 2], :);
colormap([2 4], :) = colormap([4 2], :);
colormap([3 4], :) = colormap([4 3], :);


set(pieHandle(1),  'FaceColor', colormap(1,:));
set(pieHandle(3),  'FaceColor', colormap(2,:));
set(pieHandle(5),  'FaceColor', colormap(3,:));
set(pieHandle(7),  'FaceColor', colormap(4,:));
set(pieHandle(9),  'FaceColor', colormap(5,:));
set(pieHandle(11), 'FaceColor', colormap(6,:));
set(pieHandle(13), 'FaceColor', colormap(7,:));

set(gcf, 'position', [5   278   818   653]);
%print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure3.jpg', '-r300', '-djpeg95')
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure4', '-dsvg')


% Obtain Brain ic density.
brainIdx = find(topLabelIdx==1);
brainXyz = icwiseData(brainIdx, 3:5);
brainDipoleDepth = icwiseData(brainIdx, 8);


% Obtain muscle rejection rate.
muscleIdx  = find(topLabelIdx==2);
muscleProb = icwiseData(muscleIdx,10);
muscleRejRate = sum(muscleProb>0.8)/length(muscleProb);
muscleXyz = icwiseData(muscleIdx(muscleProb>0.8), 3:5);
muscleDipoleDepth = icwiseData(muscleIdx(muscleProb>0.8), 8);


% Obtain eye rejection rate.
eyeIdx  = find(topLabelIdx==3);
eyeProb = icwiseData(eyeIdx,11);
eyeRejRate = sum(eyeProb>0.8)/length(eyeProb);
eyeXyz = icwiseData(eyeIdx(eyeProb>0.8), 3:5);
eyeDipoleDepth = icwiseData(eyeIdx(eyeProb>0.8), 8);

% Obtain heart rejection rate.
heartIdx  = find(topLabelIdx==4);
heartProb = icwiseData(heartIdx,12);
heartRejRate = sum(heartProb>0.8)/length(heartProb);
%heartXyz = icwiseData(heartIdx(heartProb>0.8), 3:5);

% Obtain line noise rejection rate.
lineNoiseIdx  = find(topLabelIdx==5);
lineNoiseProb = icwiseData(lineNoiseIdx,13);
lineNoiseRejRate = sum(lineNoiseProb>0.8)/length(lineNoiseProb);
%lineNoiseXyz = icwiseData(lineNoiseIdx(lineNoiseProb>0.8), 3:5);


% Obtain channel noise rejection rate.
channelNoiseIdx  = find(topLabelIdx==6);
channelNoiseProb = icwiseData(channelNoiseIdx,14);
channelNoiseRejRate = sum(channelNoiseProb>0.8)/length(channelNoiseProb);
% channelNoiseXyz = icwiseData(channelNoiseIdx(channelNoiseProb>0.8), 3:5);


% Obtain Cleaning Level 1 num ICs reject.
cleaningLevel1NumRejRate = (sum(muscleProb>0.8)+sum(eyeProb>0.8)+sum(heartProb>0.8)+sum(lineNoiseProb>0.8)+sum(channelNoiseProb>0.8))/...
     (sum(topLabelIdx==2)+sum(topLabelIdx==3)+sum(topLabelIdx==4)+sum(topLabelIdx==5)+sum(topLabelIdx==6));

% Plot dipole density for Muscle, Eye, and Brain.
plotDipoleDensity_axialSagittalCoronalTrio(muscleXyz, 20)
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure5.jpg', '-r300', '-djpeg95')

plotDipoleDensity_axialSagittalCoronalTrio(eyeXyz, 20)
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure6.jpg', '-r300', '-djpeg95') 

plotDipoleDensity_axialSagittalCoronalTrio(brainXyz(brainDipoleDepth<1,:), 20)
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure7.jpg', '-r300', '-djpeg95')
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison for PSD. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dummy data.
EEG = pop_loadset('filename','AAB35.set','filepath','/data/mobi/Hiroki/p3100_import_process/', 'loadmode', 'info');

freqs = EEG.Hiroki.precleanedEEG.spectraFreqs;
deltaIdx = find(freqs <= 4);
thetaIdx = find(freqs > 4 & freqs <= 8);
alphaIdx = find(freqs > 8 & freqs <= 13);
betaIdx  = find(freqs > 13 & freqs <= 30);
gammaIdx = find(freqs > 30);

psd_orig_mean   = squeeze(mean(psd_orig,1));
psd_asr_mean    = squeeze(mean(psd_asr,1));
psd_clean1_mean = squeeze(mean(psd_clean1,1));
psd_clean2_mean = squeeze(mean(psd_clean2,1));

figure
iterIdx = 0;
for cleaningIdx = 1:3
    switch cleaningIdx
%         case 1
%             currentData = psd_asr_mean-psd_orig_mean;
%         case 2
%             currentData = psd_clean1_mean-psd_orig_mean;
%         case 3
%             currentData = psd_clean2_mean-psd_orig_mean;
        case 1
            currentData = psd_orig_mean-psd_asr_mean;
        case 2
            currentData = psd_orig_mean-psd_clean1_mean;
        case 3
            currentData = psd_orig_mean-psd_clean2_mean;
    end
    
    for freqIdx = 1:5
        switch freqIdx
            case 1
                currentFreqIdx = deltaIdx;
                colorRange = [-5 5];
            case 2
                currentFreqIdx = thetaIdx;
                colorRange = [-5 5];
            case 3
                currentFreqIdx = alphaIdx;
                colorRange = [-5 5];
            case 4
                currentFreqIdx = betaIdx;
                colorRange = [-5 5];
            case 5
                currentFreqIdx = gammaIdx;
                colorRange = [-5 5];
        end
        
        iterIdx = iterIdx+1;
        subplot(4,5,iterIdx)
        
        topoplot(mean(currentData(:,currentFreqIdx),2), EEG.chanlocs, 'maplimits', colorRange)
        colorbar
    end
end

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 20)
set(gcf, 'position', [1           2        1858         929])
print('/data/mobi/Hiroki/p4100_groupLevelStructure/Figure7.jpg', '-r300', '-djpeg95')




% 
% figure
% histHandle = histogram(icwiseData(insideBrainIcIdx,8), 'FaceColor', [0.66 0.76 1]);
% 
% 
% title(sprintf('ASR rejection rate (M=%.1f, SD=%.1f, range %.0f-%.1f)',...
%               mean(1-subjectwiseCleaningStats(:,1))*100, ...
%               std(1-subjectwiseCleaningStats(:,1))*100, ...
%               min(1-subjectwiseCleaningStats(:,1))*100, ...
%               max(1-subjectwiseCleaningStats(:,1))*100))
% xlabel('Datapoint rejection rate (%)')
% ylabel('Number of datasets')
% set(findall(gcf, '-property', 'fontsize'), 'fontsize', 14)
% 
% 
% 
% 
% 
% 
% figure
% raincloud_plot(1-subjectwiseCleaningStats(:,1), 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount', 0.175, ...
%                'color', [0.66 0.76 1], 'line_width', 1, 'lwr_bnd', 0.5)
% 
% figure
% inputData = {subjectwiseCleaningStats(:,2) subjectwiseCleaningStats(:,3) subjectwiseCleaningStats(:,4)};
% rm_raincloud(inputData, [1 0 0; 0 1 0; 0 0 1], 1, 'ks')