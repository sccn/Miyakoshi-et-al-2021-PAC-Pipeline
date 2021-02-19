% 02/16/2021 Makoto. Sleep vs Awake separated. Box plot used.
% 02/08/2021 Makoto. Developed.

clear

load /data/mobi/Hiroki/p4200_select50subjects/summaryData
addpath /data/projects/makoto/Tools/RainCloudPlots-master/tutorial_matlab
addpath /data/mobi/Daisuke/code


%% 1. Histogram of ASR datapoint rejections.
figure
histHandle = histogram([1-subjectwiseCleaningStats(:,1)]*100, 0:1:40, 'FaceColor', [0.66 0.76 1]);
title(sprintf('ASR rejection rate (M=%.1f, SD=%.1f, range %.0f-%.1f)',...
              mean(1-subjectwiseCleaningStats(:,1))*100, ...
              std(1-subjectwiseCleaningStats(:,1))*100, ...
              min(1-subjectwiseCleaningStats(:,1))*100, ...
              max(1-subjectwiseCleaningStats(:,1))*100))
xlabel('Datapoint rejection rate (%)')
ylabel('Number of datasets')
ylim([0 130])
axis square
set(findall(gcf, '-property', 'fontsize'), 'fontsize', 24)
set(gcf, 'position', [1           2        1858         929])
print('/data/mobi/Hiroki/p4210_generatePlots/asrWindowRej', '-dsvg')



%% 2. Number of ICs rejected.
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
print('/data/mobi/Hiroki/p4210_generatePlots/numICsRej', '-dsvg')


%% 3. Class label rate. subjIdx, icIdx, x, y, z, singleOrDouble, rv, dipoleDepth, classLabelProbability
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
currentColormap = cbrewer('qual', 'Set3', 7);

currentColormap([2 7], :) = currentColormap([7 2], :);
currentColormap([2 6], :) = currentColormap([6 2], :);
currentColormap([2 4], :) = currentColormap([4 2], :);
currentColormap([3 4], :) = currentColormap([4 3], :);


set(pieHandle(1),  'FaceColor', currentColormap(1,:));
set(pieHandle(3),  'FaceColor', currentColormap(2,:));
set(pieHandle(5),  'FaceColor', currentColormap(3,:));
set(pieHandle(7),  'FaceColor', currentColormap(4,:));
set(pieHandle(9),  'FaceColor', currentColormap(5,:));
set(pieHandle(11), 'FaceColor', currentColormap(6,:));
set(pieHandle(13), 'FaceColor', currentColormap(7,:));

set(gcf, 'position', [5   278   818   653]);
%print('/data/mobi/Hiroki/p3200_forPoster/Figure3.jpg', '-r300', '-djpeg95')
print('/data/mobi/Hiroki/p4210_generatePlots/icLabelPieChart', '-dsvg')


%% 4. Diple depth for each class label for Level-1 IC cleaning.

    % % Obtain Brain ic density.
    % brainIdx = find(topLabelIdx==1);
    % brainXyz = icwiseData(brainIdx, 3:5);
    % brainDipoleDepth = icwiseData(brainIdx, 8);
    % 
    % % Obtain muscle rejection rate.
    % muscleIdx  = find(topLabelIdx==2);
    % muscleProb = icwiseData(muscleIdx,10);
    % muscleRejRate = sum(muscleProb>0.8)/length(muscleProb);
    % muscleXyz = icwiseData(muscleIdx(muscleProb>0.8), 3:5);
    % muscleDipoleDepth = icwiseData(muscleIdx(muscleProb>0.8), 8);
    % 
    % % Obtain eye rejection rate.
    % eyeIdx  = find(topLabelIdx==3);
    % eyeProb = icwiseData(eyeIdx,11);
    % eyeRejRate = sum(eyeProb>0.8)/length(eyeProb);
    % eyeXyz = icwiseData(eyeIdx(eyeProb>0.8), 3:5);
    % eyeDipoleDepth = icwiseData(eyeIdx(eyeProb>0.8), 8);
    % 
    % mean(brainDipoleDepth)
    % mean(muscleDipoleDepth)
    % mean(eyeDipoleDepth)
    % 
    % trimmedBrainDipDepth  = brainDipoleDepth(brainDipoleDepth>prctile(brainDipoleDepth,1) & brainDipoleDepth<prctile(brainDipoleDepth,99));
    % trimmedMuscleDipDepth = muscleDipoleDepth(muscleDipoleDepth>prctile(muscleDipoleDepth,1) & muscleDipoleDepth<prctile(muscleDipoleDepth,99));
    % trimmedEyeDipDepth    = eyeDipoleDepth(eyeDipoleDepth>prctile(eyeDipoleDepth,1) & eyeDipoleDepth<prctile(eyeDipoleDepth,99));
    % 
    % figure
    % subplot(1,3,1)
    % hist(trimmedBrainDipDepth, 50)
    % subplot(1,3,2)
    % hist(trimmedMuscleDipDepth, 50)
    % subplot(1,3,3)
    % hist(trimmedEyeDipDepth, 50)
    % 
    % 
    % 
    % 
    % % Obtain heart rejection rate.
    % heartIdx  = find(topLabelIdx==4);
    % heartProb = icwiseData(heartIdx,12);
    % heartRejRate = sum(heartProb>0.8)/length(heartProb);
    % %heartXyz = icwiseData(heartIdx(heartProb>0.8), 3:5);
    % 
    % % Obtain line noise rejection rate.
    % lineNoiseIdx  = find(topLabelIdx==5);
    % lineNoiseProb = icwiseData(lineNoiseIdx,13);
    % lineNoiseRejRate = sum(lineNoiseProb>0.8)/length(lineNoiseProb);
    % %lineNoiseXyz = icwiseData(lineNoiseIdx(lineNoiseProb>0.8), 3:5);
    % 
    % 
    % % Obtain channel noise rejection rate.
    % channelNoiseIdx  = find(topLabelIdx==6);
    % channelNoiseProb = icwiseData(channelNoiseIdx,14);
    % channelNoiseRejRate = sum(channelNoiseProb>0.8)/length(channelNoiseProb);
    % % channelNoiseXyz = icwiseData(channelNoiseIdx(channelNoiseProb>0.8), 3:5);
    % 
    % % Obtain Cleaning Level 1 num ICs reject.
    % cleaningLevel1NumRejRate = (sum(muscleProb>0.8)+sum(eyeProb>0.8)+sum(heartProb>0.8)+sum(lineNoiseProb>0.8)+sum(channelNoiseProb>0.8))/...
    %      (sum(topLabelIdx==2)+sum(topLabelIdx==3)+sum(topLabelIdx==4)+sum(topLabelIdx==5)+sum(topLabelIdx==6));

 
 
%% 5. Plot dipole density for Muscle, Eye, and Brain.

% Obtain Brain ic density.
brainIdx = find(topLabelIdx==1);
brainXyz = icwiseData(brainIdx, 3:5);
brainDipoleDepth = icwiseData(brainIdx, 8);

% Obtain muscle rejection rate.
muscleIdx = find(topLabelIdx==2);
muscleXyz = icwiseData(muscleIdx, 3:5);

% Obtain eye rejection rate.
eyeIdx = find(topLabelIdx==3);
eyeXyz = icwiseData(eyeIdx, 3:5);

% Obtain heart rejection rate.
heartIdx = find(topLabelIdx==4);
heartXyz = icwiseData(heartIdx, 3:5);

plotDipoleDensity_axialSagittalCoronalTrio(brainXyz(brainDipoleDepth<1,:), 20)
colormap('jet')
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4210_generatePlots/dipDensityBrain.jpg', '-r300', '-djpeg98')

plotDipoleDensity_axialSagittalCoronalTrio(muscleXyz, 20)
colormap('jet')
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4210_generatePlots/dipDensityMuscle.jpg', '-r300', '-djpeg98')

plotDipoleDensity_axialSagittalCoronalTrio(eyeXyz, 20)
colormap('jet')
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4210_generatePlots/dipDensityEye.jpg', '-r300', '-djpeg98') 
 
plotDipoleDensity_axialSagittalCoronalTrio(heartXyz, 20)
colormap('jet')
set(gcf, 'position', [5   278   818   653], 'invertHardcopy', 'off')
print('/data/mobi/Hiroki/p4210_generatePlots/dipDensityHeart.jpg', '-r300', '-djpeg98') 

colorbar
print('/data/mobi/Hiroki/p4210_generatePlots/colorbar.jpg', '-r300', '-djpeg98')

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Comparison for PSD. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dummy data.
EEG = pop_loadset('filename','AAB35.set','filepath','/data/mobi/Hiroki/p4000_import_preprocess_n600/', 'loadmode', 'info');

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
colorRange = [-4 4];
for cleaningIdx = 1:3
    switch cleaningIdx
        case 1
            currentData = psd_asr_mean-psd_orig_mean;
        case 2
            currentData = psd_clean1_mean-psd_orig_mean;
        case 3
            currentData = psd_clean2_mean-psd_orig_mean;
    end
    
    for freqIdx = 1:5
        switch freqIdx
            case 1
                currentFreqIdx = deltaIdx;
            case 2
                currentFreqIdx = thetaIdx;
            case 3
                currentFreqIdx = alphaIdx;
            case 4
                currentFreqIdx = betaIdx;
            case 5
                currentFreqIdx = gammaIdx;
        end
        
        iterIdx = iterIdx+1;
        subplot(3,5,iterIdx)
        
        topoplot(mean(currentData(:,currentFreqIdx),2), EEG.chanlocs, 'maplimits', colorRange)
        
        if iterIdx == 5
            currentPosition = get(gca, 'position')
            colorbar
            set(gca, 'position', currentPosition)
        end
    end
end

set(findall(gcf, '-property', 'fontsize'), 'fontsize', 20)
set(gcf, 'position', [1           2        1858         929])
print('/data/mobi/Hiroki/p4210_generatePlots/psdDiffForFreqBands.jpg', '-r300', '-djpeg98')

% figure
% raincloud_plot(1-subjectwiseCleaningStats(:,1), 'box_on', 1, 'box_dodge', 1, 'box_dodge_amount', 0.175, ...
%                'color', [0.66 0.76 1], 'line_width', 1, 'lwr_bnd', 0.5)
% 
% figure
% inputData = {subjectwiseCleaningStats(:,2) subjectwiseCleaningStats(:,3) subjectwiseCleaningStats(:,4)};
% rm_raincloud(inputData, [1 0 0; 0 1 0; 0 0 1], 1, 'ks')

%% Generate MI plots with bar graphs (original)

% Load dummy data.
EEG = pop_loadset('filename','AAB35.set','filepath','/data/mobi/Hiroki/p4000_import_preprocess_n600/', 'loadmode', 'info');
chanLabels = {EEG.chanlocs.labels}';

% Determine trimming threshold.
trimmingProb = 0.05;

% Loop for awake and sleep.
for sleepAwakeIdx = 1:2
    
    switch sleepAwakeIdx
        case 1
            currentSubjIdx  = find(pacGroupIdx==1 | pacGroupIdx==3);
            currentSubjData = pacTensor(currentSubjIdx,:,:);
        case 2
            currentSubjIdx  = find(pacGroupIdx==2 | pacGroupIdx==4);
            currentSubjData = pacTensor(currentSubjIdx,:,:);
    end
    
    sortedData  = sort(currentSubjData,1,'descend');
    trimmedData = sortedData(round(size(sortedData,1)*trimmingProb)+1:end,:,:);
    
    figure
    subplot(2,3,1)
    topoplot(squeeze(mean(zscore(trimmedData(:,:,1),0,2))), EEG.chanlocs, 'conv', 'off', 'maplimits', [-3 3])
    subplot(2,3,2)
    topoplot(squeeze(mean(zscore(trimmedData(:,:,2),0,2))), EEG.chanlocs, 'conv', 'off', 'maplimits', [-3 3])
    subplot(2,3,3)
    topoplot(squeeze(mean(zscore(trimmedData(:,:,3),0,2))), EEG.chanlocs, 'conv', 'off', 'maplimits', [-3 3])
    
    subplot(2,3,4)
    bar(mean(trimmedData(:,:,1),1), 'FaceColor', [0.66 0.76 1]);
    set(gca, 'xtick', [1:21], 'xticklabel', chanLabels)
    hold on
    errorbar(1:21, mean(trimmedData(:,:,1),1), zeros(size(trimmedData(:,:,1),2),1), std(trimmedData(:,:,1),0,1), 'color', [0 0 0], 'linestyle', 'none')
    ylabel('Canolty''s MI')
    title('Without cleaning, Mean+/-SD')
    switch sleepAwakeIdx
        case 1
            ylim([0 10])
        case 2
            ylim([0 4.5])
    end
    
    subplot(2,3,5)
    bar(mean(trimmedData(:,:,2),1), 'FaceColor', [0.66 0.76 1]);
    set(gca, 'xtick', [1:21], 'xticklabel', chanLabels)
    hold on
    errorbar(1:21, mean(trimmedData(:,:,2),1), zeros(size(trimmedData(:,:,2),2),1), std(trimmedData(:,:,2),0,1), 'color', [0 0 0], 'linestyle', 'none')
    ylabel('Canolty''s MI')
    title('ICrej 1, Mean+/-SD')
    switch sleepAwakeIdx
        case 1
            ylim([0 10])
        case 2
            ylim([0 4.5])
    end
    
    subplot(2,3,6)
    bar(mean(trimmedData(:,:,3),1), 'FaceColor', [0.66 0.76 1]);
    set(gca, 'xtick', [1:21], 'xticklabel', chanLabels)
    hold on
    errorbar(1:21, mean(trimmedData(:,:,3),1), zeros(size(trimmedData(:,:,3),2),1), std(trimmedData(:,:,3),0,1), 'color', [0 0 0], 'linestyle', 'none')
    ylabel('Canolty''s MI')
    title('ICrej 2, Mean+/-SD')
    switch sleepAwakeIdx
        case 1
            ylim([0 10])
        case 2
            ylim([0 4.5])
    end
    
    % Print.
    set(gcf, 'position', [1           2        1858         929])
    print(sprintf('/data/mobi/Hiroki/p4210_generatePlots/miTopo_%d', sleepAwakeIdx), '-dsvg')
end



%% Generate MI plots with box plots (median+quartiles, 02/16/2021)

clear

% Load summary data.
load /data/mobi/Hiroki/p4200_select50subjects/summaryData

% Load dummy data.
EEG = pop_loadset('filename','AAB35.set','filepath','/data/mobi/Hiroki/p4000_import_preprocess_n600/', 'loadmode', 'info');
chanLabels = {EEG.chanlocs.labels}';

% Determine trimming threshold. pacTensor(:,:,1:3) are MI, pacTensor(:,:,4:6) are RVL.
trimmingProb = 0.05;

% Loop for awake and sleep.
for sleepAwakeIdx = 1:2
    
    % Average Time 1 and Time 2.
    switch sleepAwakeIdx
        case 1
            currentSubjIdx1 = find(pacGroupIdx==1);
            currentSubjIdx3 = find(pacGroupIdx==3);
            currentSubjData = (pacTensor(currentSubjIdx1,:,:) + pacTensor(currentSubjIdx3,:,:))/2;
        case 2
            currentSubjIdx2 = find(pacGroupIdx==2);
            currentSubjIdx4 = find(pacGroupIdx==4);
            currentSubjData = (pacTensor(currentSubjIdx2,:,:) + pacTensor(currentSubjIdx4,:,:))/2;
    end
    
    sortedData  = sort(currentSubjData,1,'descend');
    trimmedData = sortedData(round(size(sortedData,1)*trimmingProb)+1:end,:,:);
  
    figure
    subplot(2,3,1)
    topoplot(squeeze(mean(zscore(trimmedData(:,:,1),0,2))), EEG.chanlocs, 'conv', 'off', 'maplimits', [-3 3])
    subplot(2,3,2)
    topoplot(squeeze(mean(zscore(trimmedData(:,:,2),0,2))), EEG.chanlocs, 'conv', 'off', 'maplimits', [-3 3])
    subplot(2,3,3)
    topoplot(squeeze(mean(zscore(trimmedData(:,:,3),0,2))), EEG.chanlocs, 'conv', 'off', 'maplimits', [-3 3])
    
    subplot(2,3,4)
    boxplot(currentSubjData(:,:,1), 'datalim', [0 16], 'width', 1, 'jitter', 0.25, 'plotstyle', 'compact')
    outlierTags = findobj(gca, 'tag', 'Outliers');
    set(outlierTags, 'Marker', '.', 'markerEdgeColor', [0.5 0.5 1], 'markersize', 8)    
    ylabel('Canolty''s MI')
    title('Cleaning Level 0, Mean+/-SD')
    set(gca, 'xtick', 1:21, 'xticklabel', chanLabels)
    ylim([0 16.5])

    subplot(2,3,5)
    boxplot(currentSubjData(:,:,2), 'datalim', [0 16], 'width', 1, 'jitter', 0.25, 'plotstyle', 'compact')
    outlierTags = findobj(gca, 'tag', 'Outliers');
    set(outlierTags, 'Marker', '.', 'markerEdgeColor', [0.5 0.5 1], 'markersize', 8)
    ylabel('Canolty''s MI')
    title('Cleaning Level 1, Mean+/-SD')
    set(gca, 'xtick', 1:21, 'xticklabel', chanLabels)
    ylim([0 16.5])
    
    subplot(2,3,6)
    boxplot(currentSubjData(:,:,3), 'datalim', [0 16], 'width', 1, 'jitter', 0.25, 'plotstyle', 'compact')
    outlierTags = findobj(gca, 'tag', 'Outliers');
    set(outlierTags, 'Marker', '.', 'markerEdgeColor', [0.5 0.5 1], 'markersize', 8)
    ylabel('Canolty''s MI')
    title('Cleaning Level 2, Mean+/-SD')
    set(gca, 'xtick', 1:21, 'xticklabel', chanLabels)
    ylim([0 16.5])

    % Print.
    set(gcf, 'position', [1           2        1858         929])
    print(sprintf('/data/mobi/Hiroki/p4210_generatePlots/miTopoBoxPlot_%d', sleepAwakeIdx), '-dsvg')
end



%%

% MI.
% inputData = {squeeze(miRvlTensor(1,1,:)) ...
%              squeeze(miRvlTensor(3,1,:))};           
% inputData = {squeeze(miRvlTensor(1,3,:)) ...
%              squeeze(miRvlTensor(3,3,:))};

% inputData = {squeeze(miRvlTensor(2,1,:)) ...
%              squeeze(miRvlTensor(4,1,:))};
inputData = {squeeze(miRvlTensor(2,3,:)) ...
             squeeze(miRvlTensor(4,3,:))}; 
         
% Apply trimming.
for loopIdx = 1:length(inputData)
    currentData = inputData{loopIdx};
    trimmedData = currentData(currentData<prctile(currentData,90));
    inputData{loopIdx} = trimmedData;
end
         
% Obtain fancy colors from cbrewer.
addpath('/data/projects/makoto/Tools/cbrewer/cbrewer')
currentColorMap = cbrewer('qual', 'Set3', 10);
    % tmp = 1:10;
    % figure; imagesc(tmp)
    % colormap(currentColormap)
currentColorMap = currentColorMap(4:5,:);
         
figure
rm_raincloud(inputData, currentColorMap, 1, 'ks')

figure
inputData = {squeeze(miRvlTensor(1,1,:)) ...
             squeeze(miRvlTensor(2,1,:)) ...
             squeeze(miRvlTensor(3,1,:)) ...
             squeeze(miRvlTensor(4,1,:)) ...
             squeeze(miRvlTensor(1,2,:)) ...
             squeeze(miRvlTensor(2,2,:)) ...
             squeeze(miRvlTensor(3,2,:)) ...
             squeeze(miRvlTensor(4,2,:)) ...
             squeeze(miRvlTensor(1,3,:)) ...
             squeeze(miRvlTensor(2,3,:)) ...
             squeeze(miRvlTensor(3,3,:)) ...
             squeeze(miRvlTensor(4,3,:))};  
         
% Obtain fancy colors from cbrewer.
addpath('/data/projects/makoto/Tools/cbrewer/cbrewer')
currentColormap = cbrewer('qual', 'Set3', 12);

%% Calculate state-separate asrVarRej and numIcRej 

load /data/mobi/Hiroki/p4200_select50subjects/summaryData


currentSubjIdx  = find(pacGroupIdx==1 | pacGroupIdx==3);
awakeData = subjectwiseCleaningStats(currentSubjIdx,:);

currentSubjIdx  = find(pacGroupIdx==2 | pacGroupIdx==4);
sleepData = subjectwiseCleaningStats(currentSubjIdx,:);

1-mean(awakeData); % 0.052 (0.060) 6.3 (2.6) 11.4 (2.8)
std(awakeData)     % range 0-13, 2-18.

min(awakeData(:,3:4))
max(awakeData(:,3:4))

% Awake
1-mean(awakeData); % 0.052 (0.060) 6.3 (2.6) 11.4 (2.8)
std(awakeData)     % range 0-0.33, 0-13, 2-18.

min(awakeData(:,[1 3:4]))
max(awakeData(:,[1 3:4]))

% Sleep
21-mean(sleepData); % 0.056 (0.085) 4.2 (2.4) 7.3 (3.3)
std(sleepData)      % range 0-0.39, 0-11, 0-15.

min(sleepData(:,[1 3:4]))
max(sleepData(:,[1 3:4]))

[H,P,CI,STATS] = ttest(awakeData-sleepData);

% ASR: 5.2(6.0)  0-33 vs. 5.6(8.5) 0-39  t(199)<1 n.s.
% ic1: 6.3(2.6)  0-13 vs. 4.2(2.4) 0-11, t(199)=8.6 p<0.001
% ic2: 11.4(2.8) 2-18 vs. 7.3(3.3) 0-15, t(199)=13.4 p<0.001



    
%     % 2. Obtain datapoint rejection rate by ASR.
%     subjectwiseCleaningStats(setIdx,1) = sum(EEG.etc.clean_sample_mask)/length(EEG.etc.clean_sample_mask);
%     
%     % 3. Obtain the number of ICs left for three levels of cleaning.
%     subjectwiseCleaningStats(setIdx,2) = size(EEG.icaweights,1);
%     subjectwiseCleaningStats(setIdx,3) = size(EEG.Hiroki.icRej1.icaweights,1);
%     subjectwiseCleaningStats(setIdx,4) = size(EEG.Hiroki.icRej2.icaweights,1);