% 01/21/2020 Makoto. Debugged.
% 11/04/2020 Makoto. ASR with window rejection. PAC with peri-edge window rejection. 
% 02/28/2020 Makoto. Notch filter supported.
% 02/13/2020 Makoto. Backward compatibility to pre-cleaned data supported.
% 02/12/2020 Makoto. Channel process updated.
% 02/07/2020 Makoto. Written for Hiroki.

currentTic = tic;

eeglabPath = which('eeglab');
eeglabPath = eeglabPath(1:end-9);
addpath([eeglabPath '/plugins/Fieldtrip-lite20200130/fileio']); % Only for my environement. Bug?
allEdf = dir('/data/mobi/Hiroki/p3000_raw/*.edf');
for edfIdx = 1 %:length(allEdf)
    
    try
        % Import the .edf file.
        currentEdfName = allEdf(edfIdx).name;
        EEG = pop_biosig(['/data/mobi/Hiroki/p3000_raw/' currentEdfName], 'blockrange',[60 360]);
        EEG.setname = currentEdfName(1:end-4);
        EEG.subject = currentEdfName(1:end-4);
        
        % If sampling rate is not 200, downsample it to 200 Hz.
        if EEG.srate ~= 200
            EEG = pop_resample(EEG, 200);
        end
        
%         % Trim data to 1-6 min.
%         if EEG.xmax > 360
%             EEG = pop_select(EEG, 'time', [60 360]);
%         end
        
        % Reject non-EEG channels.
        targetChannelList = {'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F7' 'F8' 'T3' 'T4' 'T5' 'T6' 'FZ' 'CZ' 'PZ' 'A1' 'A2'};
        channelLabels = {EEG.chanlocs.labels}';
        identifiedChList = zeros(length(targetChannelList),1);
        for chIdx = 1:length(targetChannelList)
            currentIdx = find(contains(channelLabels, targetChannelList(chIdx), 'ignorecase', true));
            identifiedChList(chIdx) = currentIdx(1);
        end
        EEG = pop_select(EEG, 'channel', identifiedChList);
        
        % Rename channels.
        selectedChannelLabels = {EEG.chanlocs.labels}';
        sortingIdxList = zeros(length(targetChannelList),1);
        for chIdx = 1:length(targetChannelList)
            currentIdx = find(contains(selectedChannelLabels, targetChannelList(chIdx), 'ignorecase', true));
            sortingIdxList(chIdx) = currentIdx;
            EEG.chanlocs(currentIdx).labels = targetChannelList{chIdx};
        end
        
        % Sort channels.
        EEG.data = EEG.data(sortingIdxList,:);
        EEG.chanlocs = EEG.chanlocs(sortingIdxList);
        
        % Import channel locations.
        EEG = pop_chanedit(EEG, 'lookup', '/data/projects/makoto/Tools/eeglab14_1_2b/plugins/dipfit3.3/standard_BEM/elec/standard_1005.elc');
        chanCoordinates = [[EEG.chanlocs.X]' [EEG.chanlocs.Y]' [EEG.chanlocs.Z]'];
        if isempty(chanCoordinates)
            error('Failed to import template channel locations.');
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Compute PACT on raw data. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform high-pass filter (passband edge: 0.5 Hz)
        EEG = pop_eegfiltnew(EEG, [], 0.5, 1320, 1, [], 0);
        
        % Step 6: Run cleanLineNoise (requires PREP Pipeline toolbox) The original CleanLine was replaced (04/25/2020).
        signal      = struct('data', EEG.data, 'srate', EEG.srate);
        lineNoiseIn = struct('lineNoiseMethod', 'clean', ...
                        'lineNoiseChannels', 1:EEG.nbchan,...
                        'Fs', EEG.srate, ...
                        'lineFrequencies', 60,...
                        'p', 0.01, ...
                        'fScanBandWidth', 2, ...
                        'taperBandWidth', 2, ...
                        'taperWindowSize', 4, ...
                        'taperWindowStep', 1, ...
                        'tau', 100, ...
                        'pad', 2, ...
                        'fPassBand', [0 EEG.srate/2], ...
                        'maximumIterations', 10);
        [clnOutput, lineNoiseOut] = cleanLineNoise(signal, lineNoiseIn);
        EEG.data = clnOutput.data;
        
%         % Perform notch filter (passband edge: 0.5 Hz)
%         EEG = pop_firws(EEG, 'fcutoff', [59 61], 'ftype', 'bandstop', 'wtype', 'hamming', 'forder', 660, 'minphase', 0);
        
        % PACT on non-cleaned data.
        EEG = pac_man(EEG, 'lfoPhase', [3  4], 'hfoAmp', [35  70], 'hfoTopRatio', 2, 'hfoPool', 1, 'whichMarker', 1, 'windowLength', [], 'alpha', [1.000000e-02], 'numSurro', 2000, 'numPhaseBin', 36);
              
        % Compute PSD of the raw data.
        [spectra,freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);

        % Store the raw data.
        EEG.Hiroki.precleanedEEG.data         = single(EEG.data);
        EEG.Hiroki.precleanedEEG.spectra      = single(spectra);
        EEG.Hiroki.precleanedEEG.spectraFreqs = single(freqs);
        EEG.Hiroki.precleanedEEG.pac          = EEG.pac;
        EEG = rmfield(EEG, 'pac');

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Compute ASR and ICA. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ASR without changing data length.
%         EEG = clean_rawdata(EEG, 'off', 'off', 'off', 'off', 20, 'off', 'availableRAM_GB', 8);
        EEG = clean_rawdata(EEG, 'off', 'off', 'off', 'off', 20, 0.25, 'availableRAM_GB', 8);
        
        % ICA.
        EEG_forICA = pop_resample(EEG, 100);
        EEG_forICA = pop_runica(EEG_forICA, 'extended',1,'interupt','off');
        EEG.icaweights = EEG_forICA.icaweights;
        EEG.icasphere  = EEG_forICA.icasphere;
        EEG = eeg_checkset(EEG, 'ica');
        
        % Dipfit.
        EEG = pop_dipfit_settings( EEG, 'hdmfile', '/data/projects/makoto/Tools/eeglab14_1_2b/plugins/dipfit3.3/standard_BEM/standard_vol.mat',...
            'mrifile', '/data/projects/makoto/Tools/eeglab14_1_2b/plugins/dipfit3.3/standard_BEM/standard_mri.mat',...
            'chanfile','/data/projects/makoto/Tools/eeglab14_1_2b/plugins/dipfit3.3/standard_BEM/elec/standard_1005.elc',...
            'coord_transform', [0 0 0 0 0 -1.5708 1 1 1] , 'chansel', 1:EEG.nbchan, 'coordformat','MNI');
        EEG = pop_multifit(EEG, [1:21] ,'threshold',100,'plotopt',{'normlen' 'on'});
        EEG = fitTwoDipoles(EEG, 'LRR', 35);
        
        % ICLabel.
        EEG = pop_iclabel(EEG, 'default');
        
        % Compute PSD of the ASD-ed data.
        [spectra,freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);

        % Store ICA and dipfit results.
        EEG.Hiroki.asrCleaned.spectra      = single(spectra);
        EEG.Hiroki.asrCleaned.spectraFreqs = single(freqs);
%         EEG.Hiroki.asrCleaned.data         = single(EEG.data);
%         EEG.Hiroki.asrCleaned.icaweights   = EEG.icaweights; 
%         EEG.Hiroki.asrCleaned.icasphere    = EEG.icasphere;
%         EEG.Hiroki.asrCleaned.dipfit       = dipfit;
        
        % keep the Original ICA-ed data.
        EEG_originalICA = EEG;
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Gentle rejection. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform IC rejection using ICLabel scores and r.v. from dipole fitting.
        muscleIdx = find(EEG.etc.ic_classification.ICLabel.classifications(:,2) >= 0.80);
        eyeIdx    = find(EEG.etc.ic_classification.ICLabel.classifications(:,3) >= 0.80);
        heartIdx  = find(EEG.etc.ic_classification.ICLabel.classifications(:,4) >= 0.80);
        allArtifactIdx = unique([muscleIdx; eyeIdx; heartIdx]);
        brainIdx  = setdiff(1:EEG.nbchan, allArtifactIdx);
        rvList    = [EEG.dipfit.model.rv];
        goodRvIdx = find(rvList < 0.25); % < 20% residual variance == good ICs.
        goodIcIdx = intersect(brainIdx, goodRvIdx);
        EEG = pop_subcomp(EEG, goodIcIdx, 0, 1);
        %EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodIcIdx,:);
        
        % PACT on cleaned data.
        EEG = pac_man(EEG, 'lfoPhase', [3  4], 'hfoAmp', [35  70], 'hfoTopRatio', 2, 'hfoPool', 1, 'whichMarker', 1, 'windowLength', [], 'alpha', [1.000000e-02], 'numSurro', 2000, 'numPhaseBin', 36);

        % Compute PSD of the ASD+ICrej1 data.
        [spectra,freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);
        
        % Store ICA and dipfit results.
        EEG.Hiroki.icRej1.data         = single(EEG.data);
        EEG.Hiroki.icRej1.spectra      = single(spectra);
        EEG.Hiroki.icRej1.spectraFreqs = single(freqs);
        EEG.Hiroki.icRej1.icaweights   = EEG.icaweights; 
        EEG.Hiroki.icRej1.icasphere    = EEG.icasphere;
        EEG.Hiroki.icRej1.dipfit       = EEG.dipfit;
        EEG.Hiroki.icRej1.pac          = EEG.pac;
        EEG = rmfield(EEG, 'pac');
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Aggressive rejection. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Recover the origianl ICA-ed data.
        EEG.data       = EEG_originalICA.data;
        EEG.icaweights = EEG_originalICA.icaweights;
        EEG.icasphere  = EEG_originalICA.icasphere;
        EEG.dipfit     = EEG_originalICA.dipfit;
        EEG.icaact     = [];
        EEG.icawinv    = [];
        EEG = eeg_checkset(EEG, 'ica');
        
        % Perform IC rejection using ICLabel scores and r.v. from dipole fitting.
        % brainIdx = find(EEG.etc.ic_classification.ICLabel.classifications(:,1) >= 0.70);
        %
        % Option 2. Use Brain criterion regardless of the label probability.
        icLabelMatrix  = EEG.etc.ic_classification.ICLabel.classifications;
        icLabelClasses = EEG.etc.ic_classification.ICLabel.classes;
        brainIdx = find(strcmp(icLabelClasses, 'Brain'));
        [~,labelVector] = max(icLabelMatrix, [], 2);
        brainIdx = find(labelVector==brainIdx);
        
        rvList    = [EEG.dipfit.model.rv];
        goodRvIdx = find(rvList < 0.25); % < 15% residual variance == good ICs.
        goodIcIdx = intersect(brainIdx, goodRvIdx);
        EEG = pop_subcomp(EEG, goodIcIdx, 0, 1);
        %EEG.etc.ic_classification.ICLabel.classifications = EEG.etc.ic_classification.ICLabel.classifications(goodIcIdx,:);
        
        % PACT on cleaned data.
        EEG = pac_man(EEG, 'lfoPhase', [3  4], 'hfoAmp', [35  70], 'hfoTopRatio', 2, 'hfoPool', 1, 'whichMarker', 1, 'windowLength', [], 'alpha', [1.000000e-02], 'numSurro', 2000, 'numPhaseBin', 36);

        % Compute PSD of the ASD+ICrej1 data.
        [spectra,freqs] = spectopo(EEG.data, EEG.pnts, EEG.srate, 'overlap', 50, 'plot', 'off', 'freqfac', 4);
        
        % Store ICA and dipfit results.
        EEG.Hiroki.icRej2.data         = single(EEG.data);
        EEG.Hiroki.icRej2.spectra      = single(spectra);
        EEG.Hiroki.icRej2.spectraFreqs = single(freqs);
        EEG.Hiroki.icRej2.icaweights   = EEG.icaweights; 
        EEG.Hiroki.icRej2.icasphere    = EEG.icasphere;
        EEG.Hiroki.icRej2.dipfit       = EEG.dipfit;
        EEG.Hiroki.icRej2.pac          = EEG.pac;
        EEG = rmfield(EEG, 'pac');  
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Recover no-IC-rejected data for saving. %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Recover the origianl ICA-ed data.
        EEG.data       = EEG_originalICA.data;
        EEG.icaweights = EEG_originalICA.icaweights;
        EEG.icasphere  = EEG_originalICA.icasphere;
        EEG.dipfit     = EEG_originalICA.dipfit;
        EEG.icaact     = [];
        EEG.icawinv    = [];
        %EEG.etc.ic_classification.ICLabel = EEG_originalICA.etc.ic_classification.ICLabel;
        EEG = eeg_checkset(EEG, 'ica');
        
        % Save data.
        pop_saveset(EEG, 'filename', currentEdfName(1:end-4), 'filepath', '/data/mobi/Hiroki/p4000_import_preprocess_n600');
        
    catch
        dummyData = 0;
        save(['/data/mobi/Hiroki/p4000_import_preprocess_n600/FAILED_' currentEdfName(1:end-4) '.mat'], 'dummyData');
    end
end

toc(currentTic)

%% Validate channel names and channel order.
        % allSets = dir('/data/mobi/Hiroki/p0100_import_process/*.set');
        % allChLabels = cell(21, length(allSets));
        % for setIdx = 1:length(allSets)
        %     
        %     % Import the .edf file.
        %     loadName = allSets(setIdx).name;
        %     EEG = pop_loadset('filename', loadName, 'filepath', '/data/mobi/Hiroki/p0100_import_process', 'loadmode', 'info');
        %     
        %     % Show the channel label table.
        %     allChLabels(:,setIdx) = {EEG.chanlocs.labels}';
        % end    