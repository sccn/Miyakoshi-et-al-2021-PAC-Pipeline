% 01/19/2021 Makoto. Created.

load /data/mobi/Hiroki/p4100_groupLevelStructure/groupLevelStructure.mat
addpath /data/projects/makoto/Tools/RainCloudPlots-master/tutorial_matlab
addpath /data/mobi/Daisuke/code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute freq-band separated PSD. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load dummy data.
EEG = pop_loadset('filename','AAB35.set','filepath','/data/mobi/Hiroki/p4000_import_preprocess_n600/', 'loadmode', 'info');

freqs = EEG.Hiroki.precleanedEEG.spectraFreqs;
deltaIdx = find(freqs <= 4);
thetaIdx = find(freqs > 4 & freqs <= 8);
alphaIdx = find(freqs > 8 & freqs <= 13);
betaIdx  = find(freqs > 13 & freqs <= 30);
gammaIdx = find(freqs > 30);

origDelta = mean(psd_orig(  :,:,deltaIdx),3);
asrDelta  = mean(psd_asr(   :,:,deltaIdx),3);
cl1Delta  = mean(psd_clean1(:,:,deltaIdx),3);
cl2Delta  = mean(psd_clean2(:,:,deltaIdx),3);

origTheta = mean(psd_orig(  :,:,thetaIdx),3);
asrTheta  = mean(psd_asr(   :,:,thetaIdx),3);
cl1Theta  = mean(psd_clean1(:,:,thetaIdx),3);
cl2Theta  = mean(psd_clean2(:,:,thetaIdx),3);

origAlpha = mean(psd_orig(  :,:,alphaIdx),3);
asrAlpha  = mean(psd_asr(   :,:,alphaIdx),3);
cl1Alpha  = mean(psd_clean1(:,:,alphaIdx),3);
cl2Alpha  = mean(psd_clean2(:,:,alphaIdx),3);

origBeta = mean(psd_orig(  :,:,betaIdx),3);
asrBeta  = mean(psd_asr(   :,:,betaIdx),3);
cl1Beta  = mean(psd_clean1(:,:,betaIdx),3);
cl2Beta  = mean(psd_clean2(:,:,betaIdx),3);

origGamma = mean(psd_orig(  :,:,gammaIdx),3);
asrGamma  = mean(psd_asr(   :,:,gammaIdx),3);
cl1Gamma  = mean(psd_clean1(:,:,gammaIdx),3);
cl2Gamma  = mean(psd_clean2(:,:,gammaIdx),3);

asrDelta = origDelta-asrDelta;
cl1Delta = origDelta-cl1Delta;
cl2Delta = origDelta-cl2Delta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Before-After comparison: ASR. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,5,1)
imagesc(asrDelta<0)
title('Delta'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,2)
imagesc(asrTheta<0)
title('Theta'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,3)
imagesc(asrAlpha<0)
title('Alpha'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,4)
imagesc(asrBeta<0)
title('Beta');  xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,5)
imagesc(asrGamma<0)
title('Gamma'); xlabel('Electrode index'); ylabel('Dataset index')

suptitle('ASR Before-After subtraction (Blue, >0; Yellow, <0)')
set(gcf, 'position', [1 1 1858 929])
set(findall(gcf,'-property', 'fontsize'),'fontsize', 14)
print /data/mobi/Hiroki/p4110_powerIncreaseTest/beforeAfterAsr -djpeg98 -r72


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Before-After comparison: ASR+ICArej1. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,5,1)
imagesc(cl1Delta<0)
title('Delta'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,2)
imagesc(cl1Theta<0)
title('Theta'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,3)
imagesc(cl1Alpha<0)
title('Alpha'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,4)
imagesc(cl1Beta<0)
title('Beta');  xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,5)
imagesc(cl1Gamma<0)
title('Gamma'); xlabel('Electrode index'); ylabel('Dataset index')

suptitle('ASR+ICrej1 Before-After subtraction (Blue, >0; Yellow, <0)')
set(gcf, 'position', [1 1 1858 929])
set(findall(gcf,'-property', 'fontsize'),'fontsize', 14)
print /data/mobi/Hiroki/p4110_powerIncreaseTest/beforeAfterICrej1 -djpeg98 -r72

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Before-After comparison: ASR+ICArej2. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(1,5,1)
imagesc(cl2Delta<0)
title('Delta'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,2)
imagesc(cl2Theta<0)
title('Theta'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,3)
imagesc(cl2Alpha<0)
title('Alpha'); xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,4)
imagesc(cl2Beta<0)
title('Beta');  xlabel('Electrode index'); ylabel('Dataset index')
subplot(1,5,5)
imagesc(cl2Gamma<0)
title('Gamma'); xlabel('Electrode index'); ylabel('Dataset index')

suptitle('ASR+ICrej2 Before-After subtraction (Blue, >0; Yellow, <0)')
set(gcf, 'position', [1 1 1858 929])
set(findall(gcf,'-property', 'fontsize'),'fontsize', 14)
print /data/mobi/Hiroki/p4110_powerIncreaseTest/beforeAfterICrej2 -djpeg98 -r72


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Visualize the results. %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% subplot(1,3,1)
% imagesc(asrDelta<0)
% title('Delta ASR'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,2)
% imagesc(cl1Delta<0)
% title('Delta ASR+ICrej1'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,3)
% imagesc(cl2Delta<0)
% title('Delta ASR+ICrej2'); xlabel('Electrode index'); ylabel('Dataset index')
% 
% figure
% subplot(1,3,1)
% imagesc(asrTheta<0)
% title('Theta ASR'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,2)
% imagesc(cl1Theta<0)
% title('Theta ASR+ICrej1'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,3)
% imagesc(cl2Theta<0)
% title('Theta ASR+ICrej2'); xlabel('Electrode index'); ylabel('Dataset index')
% 
% figure
% subplot(1,3,1)
% imagesc(asrAlpha<0)
% title('Alpha ASR'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,2)
% imagesc(cl1Alpha<0)
% title('Alpha ASR+ICrej1'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,3)
% imagesc(cl2Alpha<0)
% title('Alpha ASR+ICrej2'); xlabel('Electrode index'); ylabel('Dataset index')
% 
% figure
% subplot(1,3,1)
% imagesc(asrBeta<0)
% title('Beta ASR'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,2)
% imagesc(cl1Beta<0)
% title('Beta ASR+ICrej1'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,3)
% imagesc(cl2Beta<0)
% title('Beta ASR+ICrej2'); xlabel('Electrode index'); ylabel('Dataset index')
% 
% figure
% subplot(1,3,1)
% imagesc(asrGamma<0)
% title('Gamma ASR'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,2)
% imagesc(cl1Gamma<0)
% title('Gamma ASR+ICrej1'); xlabel('Electrode index'); ylabel('Dataset index')
% subplot(1,3,3)
% imagesc(cl2Gamma<0)
% title('Gamma ASR+ICrej2'); xlabel('Electrode index'); ylabel('Dataset index')



%% Why PSD increase?

allSets = dir('/data/mobi/Hiroki/p3100_import_process/*.set')

for setIdx = 3 %1:length(allSets)
    
    loadName = allSets(setIdx).name;
    EEG = pop_loadset('filename',loadName,'filepath','/data/mobi/Hiroki/p3100_import_process/', 'loadmode', 'info');
    
    spectraTensor = cat(3, EEG.Hiroki.precleanedEEG.spectra, EEG.Hiroki.asrCleaned.spectra, EEG.Hiroki.icRej1.spectra, EEG.Hiroki.icRej2.spectra);
    specFreqs     = EEG.Hiroki.precleanedEEG.spectraFreqs;
    
    figure; set(gcf, 'position', [5   511   620   420])
    for elecIdx = 1:21
        subplot(3,7,elecIdx)
        plot(specFreqs, spectraTensor(elecIdx,:,1))
        hold on
        plot(specFreqs, spectraTensor(elecIdx,:,2))
        title(EEG.chanlocs(elecIdx).labels)
        ylim([-20 20])
    end
    suptitle('Raw-ASR')
    
    figure; set(gcf, 'position', [633   511   620   420])
    for elecIdx = 1:21
        subplot(3,7,elecIdx)
        plot(specFreqs, spectraTensor(elecIdx,:,1))
        hold on
        plot(specFreqs, spectraTensor(elecIdx,:,3))
        title(EEG.chanlocs(elecIdx).labels)
        ylim([-20 20])
    end
    suptitle('Raw-ICrej1')
    
    figure; set(gcf, 'position', [1236  511  619  420])
    for elecIdx = 1:21
        subplot(3,7,elecIdx)
        plot(specFreqs, spectraTensor(elecIdx,:,1))
        hold on
        plot(specFreqs, spectraTensor(elecIdx,:,4))
        title(EEG.chanlocs(elecIdx).labels)
        ylim([-20 20])
    end
    suptitle('Raw-ICrej2')
    
    close all
end

%% Comparebefore and after IC rejection

%allSets = dir('/data/mobi/Hiroki/p3100_import_process/*.set')
allSets = dir('/data/mobi/Hiroki/p4000_import_preprocess_n600/*.set')


for setIdx = 3 %1:length(allSets)
    
    % Load data.
    loadName = allSets(setIdx).name;
    %EEG = pop_loadset('filename',loadName,'filepath','/data/mobi/Hiroki/p3100_import_process/');
    EEG = pop_loadset('filename',loadName,'filepath','/data/mobi/Hiroki/p4000_import_preprocess_n600/');

    originalICs = EEG.icaweights;
    
    % Compute r.v. list
    rvList    = [EEG.dipfit.model.rv];
    goodRvIdx = find(rvList < 0.25); % < 20% residual variance == good ICs.
    
    
    % Determine which ICs were rejected in the Level-1 ICrej.
    rej1ICs = EEG.Hiroki.icRej1.icaweights;
    rej1SelectedIcIdx = zeros(size(rej1ICs,1),1);
    for icIdx = 1:size(rej1ICs,1)
        corrMatrix = corrcoef([originalICs; rej1ICs(icIdx,:)]');
        matchingIcIdx = find(corrMatrix(end,1:end-1)>0.99);
        rej1SelectedIcIdx(icIdx) = matchingIcIdx;
    end
    rejectedICs1 = setdiff(1:size(originalICs,1), rej1SelectedIcIdx);

    % Determine which ICs were rejected in the Level-2 ICrej.
    rej2ICs = EEG.Hiroki.icRej2.icaweights;
    rej2SelectedIcIdx = zeros(size(rej2ICs,1),1);
    for icIdx = 1:size(rej2ICs,1)
        corrMatrix = corrcoef([originalICs; rej2ICs(icIdx,:)]');
        matchingIcIdx = find(corrMatrix(end,1:end-1)>0.99);
        rej2SelectedIcIdx(icIdx) = matchingIcIdx;
    end
    rejectedICs2 = setdiff(1:size(originalICs,1), rej2SelectedIcIdx);

    % PSD evaluation.
    [spectra1,freqs1] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);

    % Perform IC rejection and PSD evaluation.
    EEG2 = pop_subcomp(EEG, rej1SelectedIcIdx, 0, 1);
    [spectra2,freqs2] = spectopo(EEG2.data, EEG2.pnts, EEG2.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);

    % PSD evaluation again.
    EEG3 = pop_subcomp(EEG, rej2SelectedIcIdx, 0, 1);
    [spectra3,freqs3] = spectopo(EEG3.data, EEG3.pnts, EEG3.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);


    % Plot 21-ch PSDs.
    figure
    for elecIdx = 1:21
        subplot(3,7,elecIdx)
        plot(EEG.Hiroki.icRej1.spectraFreqs, spectra3(elecIdx,:))
        hold on
        plot(EEG.Hiroki.icRej1.spectraFreqs, spectra2(elecIdx,:))
        plot(EEG.Hiroki.icRej1.spectraFreqs, spectra1(elecIdx,:))
        xlim([1 100])
        title(EEG.chanlocs(elecIdx).labels)
    end
    legend({'ICrej2' 'ICrej1' 'No ICrej'})

    
    % Plot 21-ch 10-s data.
    EEG_hpf  = pop_eegfiltnew(EEG,  [],20,132,1,[],0); % High-pass filter at 20 Hz.
    EEG2_hpf = pop_eegfiltnew(EEG2, [],20,132,1,[],0); % High-pass filter at 20 Hz.
    EEG3_hpf = pop_eegfiltnew(EEG3, [],20,132,1,[],0); % High-pass filter at 20 Hz.
    
    
    figure
    for elecIdx = 1:21
        subplot(3,7,elecIdx)
        p1Handle = plot(EEG3.times(20000:20101)/1000, EEG3_hpf.data(elecIdx,20000:20101));
        hold on
        p2Handle = plot(EEG2.times(20000:20101)/1000, EEG2_hpf.data(elecIdx,20000:20101));
        p3Handle = plot(EEG.times( 20000:20101)/1000, EEG_hpf.data( elecIdx,20000:20101));
        xlabel('Time (s)')
        xlim([100 100.5])
        ylim([-50 50])
        title(EEG.chanlocs(elecIdx).labels)
    end
    legend({'ICrej2' 'ICrej1' 'No ICrej'})

    
%     % Plot rejection-dependent power change.
%     for rejIcIdxIdx = 1:length(rejectedICs2)
%         
%         rejIcIdx = 
%         
%         
%         
%         
%         
%     end
end