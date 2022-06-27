%%
[y, Fs] = audioread('myvoice9.wav'); sound(y, Fs);
% figure; plot(y);
%%
Fs = 44100;
T = 10;
t = [0:1/Fs:T].';
f = [800]; %114*[1 2 3 4 5 6 7];
ys = sin(2*pi*f.*t);
y = sum(ys, 2);
sound(y./length(f), Fs);
% figure; plot(y(1:1024)./length(f))

% lowpass(y(1:2048), 2000, Fs)
% figure; plot(t, ys); hold on; plot(t, y, 'k', 'LineWidth', 2);
% figure; plot(abs(fft(y(1200+[1:1024]))));
% figure; plot(y);
%%
clear sound 
%% load data
tic
clearvars -except y; close all;
% bin_filename = '2tx_4rx_single_tone/adc_data_440HzTone_speaker0.93m.bin';
% copyfile(bin_filename, 'C:\ti\mmwave_studio_02_01_01_00\mmWaveStudio\PostProc\adc_data.bin');

rawDataReader('0305/0305.setup.json','', 'radarCube', 0)
load('radarCube.mat');

% save(sprintf('%s.mat', load_filename(1:end-4)), 'radarCube');
% verify_data('radarCube.mat');
% load('2tx_4rx_single_tone_0211/adc_data_440HzTone_door0.7m_trashbag0.93m_speaker1.58m.mat');
% load('2tx_4rx_single_tone_0211/adc_data_440HzTone_trashbag0.93m_speaker1.58m.mat');
% load('0305/multitone/radarCube_traderjoebag_nine.mat');

% transform radarCube into a better representation
params.startFreq = radarCube.rfParams.startFreq;
params.freqSlope = radarCube.rfParams.freqSlope;
params.sampleRate = radarCube.rfParams.sampleRate;
params.bandwidth = radarCube.rfParams.bandwidth;
params.rangeResolutionsInMeters = radarCube.rfParams.rangeResolutionsInMeters;
params.dopplerResolutionMps = radarCube.rfParams.dopplerResolutionMps;
params.framePeriodicity = radarCube.rfParams.framePeriodicity;
params.numFrames = radarCube.dim.numFrames;
params.numChirps = radarCube.dim.numChirps;
params.numRxChan = radarCube.dim.numRxChan;
params.numSamplePerChirp = radarCube.dim.numRangeBins;
params.opRangeFFTSize = max(2^14, 2^nextpow2(params.numSamplePerChirp));
params.opDopplerFFTSize = max(2^10, 2^nextpow2(params.numChirps));
params.chirpCycleTime = radarCube.rfParams.chirpCycleTime;

% params.w_steps = 1000;
% params.w_stepsize = 2*pi/(params.w_steps*params.opRangeFFTSize);
% params.w =[0:params.w_steps*params.opRangeFFTSize-1]*params.w_stepsize;
% params.chirpSincWindowSize = params.numSamplePerChirp;
% params.chirpSincs = exp(-1j*(params.w)*(params.chirpSincWindowSize-1)/2).*...
%                     diric(params.w, params.chirpSincWindowSize);
% params.subChirpSincWindowSize = 256;   
% params.subChirpSincs = exp(-1j*(params.w)*(params.subChirpSincWindowSize-1)/2).*...
%                     diric(params.w, params.subChirpSincWindowSize);
% params.chirpSearchSincHashMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
% params.chirpSincHashMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
% params.subChirpSincHashMap = containers.Map('KeyType', 'double', 'ValueType', 'any');
               
% w_diff = params.w_stepsize*[-params.w_steps:params.w_steps];
% params.sincs_adjust = 1./...
%     (exp(-1j*(w_diff)*(params.sincWindowSize-1)/2).*...
%     sin(params.sincWindowSize*w_diff/2)./(params.sincWindowSize*sin(w_diff/2)));
%                 
%                 
datacube.params = params;
% datacube.adcdata = zeros(params.numChirps, params.numRxChan, params.);
for ii=1:length(radarCube.data)
    frame = radarCube.data{ii};
    frame = ifft(frame, [], 3);
%     for jj = 1:params.numChirps
%         for kk = 1:params.numRxChan
%             adc_data = ifft(squeeze(frame(jj,kk,:)));
%             frame(jj,kk,:) = adc_data;
%         end
%     end
    datacube.adcdata{ii} = frame;
end
% clear radarCube
params

% figure; 
% for ii=1:params.numFrames
%     subplot(params.numFrames,1,ii);
%     plot(real(squeeze(datacube.adcdata{ii}(1,1,:)))); hold on;
%     plot(imag(squeeze(datacube.adcdata{ii}(1,1,:))));
% end
toc
%% post processing -- process only one frame
selected_frame = 1;
selected_chirp = 1;
selected_rxchain = 1; 
plot_range_fft = 1;
plot_doppler_fft = 1;
plot_phase_across_chirps = 1;
cf = 0;

% ragne FFT
if plot_range_fft
    frame = datacube.adcdata{selected_frame};
    adc_data = squeeze(frame(selected_chirp,selected_rxchain,:));
    rangeFFT(adc_data, datacube.params);
%     save('adc_data.mat', 'adc_data');
%     figure;
%     rangefft_output = fft(adc_data, 1024*8);
%     findpeaks(abs(rangefft_output), 'MinPeakProminence',3e5)
%     [pks,locs,w,p] = findpeaks(abs(rangefft_output), 'MinPeakProminence',3e5);
end

%%
close all;
peak = 142 + [-2:2];
phases = zeros(length(peak), params.numChirps);
for ii=1:params.numChirps
    adc_data = double(squeeze(datacube.adcdata{1}(ii,1,:))).';
    rangefft = fft(adc_data, params.opRangeFFTSize);
    phases(:, ii) = angle(rangefft(peak));
end
polyFitOrder = 9;
[vibration_spectrum1, corrected_phases, I] = phases_to_spectrum(phases, params, polyFitOrder, 1);
% I = round(400/(1/params.chirpCycleTime/params.opDopplerFFTSize), 0);
% abs(vibration_spectrum(3,I+1))
% noise_sample_idx = floor(params.opDopplerFFTSize/2*0.8):params.opDopplerFFTSize/2;
% db(abs(vibration_spectrum(3,I+1))/mean(abs(vibration_spectrum(3,noise_sample_idx))))


% figure; plot(corrected_phases(3,:));
% figure(55); highpass(y(1:2048), 500, Fs)
phases = corrected_phases(3,:);
figure(56); bandpass(phases,[220 1500], 1/params.chirpCycleTime, 'Steepness', 0.9)
% bpFilt = designfilt('bandpassiir','FilterOrder',10, 'HalfPowerFrequency1',245,'HalfPowerFrequency2',265,...
%     'DesignMethod','butter', 'SampleRate',14925);
% phases2 = filtfilt(bpFilt, phases);

bpFilt2 = designfilt('bandpassfir','FilterOrder',2000, 'CutoffFrequency1',225,'CutoffFrequency2',235,...
    'SampleRate',14925);
phases2 = conv(bpFilt2.Coefficients, phases); phases2 = phases2(1+1000:end-1000);
figure; subplot(2,1,1); plot(phases); hold on; plot(phases2);
subplot(2,1,2); plot(db(fft(phases, 4096))); hold on; plot(db(fft(phases2, 4096)));

fundamental = 114;
eq_phases = phases; %zeros(size(phases));
eq_db = [0 20 27 32 34 37 40];
for ii=1:6
    bpFilt2 = designfilt('bandpassfir','FilterOrder',2000, 'CutoffFrequency1',ii*fundamental-5,...
        'CutoffFrequency2',ii*fundamental+5, 'SampleRate',14925);
    phases2 = conv(bpFilt2.Coefficients, phases); phases2 = phases2(1+1000:end-1000);
    eq_phases = eq_phases + 10^(eq_db(ii)/20) * phases2 - phases2;
end
figure(58); bandpass(eq_phases,[220 235], 1/params.chirpCycleTime)
%%
peak=142;
phases = zeros(params.numFrames, params.numChirps, 3);
vibration_spectrum = zeros(params.numFrames, params.numChirps, 3);
corrected_phases = zeros(params.numFrames, params.numChirps, 3);
for kk=1:3
    for ii=1:params.numFrames
        for jj=1:params.numChirps
            adc_data = double(squeeze(datacube.adcdata{ii}(jj,kk,:))).';
            rangefft = fft(adc_data, params.opRangeFFTSize);
    %         phases((ii-1)*params.numChirps+jj) = angle(rangefft(peak));
            phases(ii, jj, kk) = angle(rangefft(peak));
        end
    end
    polyFitOrder = 9;
    [vibration_spectrum(:,:,kk), corrected_phases(:,:,kk), I] = phases_to_spectrum(phases(:,:,kk), params, polyFitOrder, 0);
end

% figure; plot(1:2048, squeeze(phases(1,:,:)) );
% figure; plot(1:2048, abs(squeeze(phases(1,:,:)) - squeeze(phases(1,:,1).') ));
% figure; plot(1:2048, squeeze(corrected_phases(12,:,:)));
% 
% figure;
% subplot(2,1,1);
% for ii=1:3
%     plot(db(vibration_spectrum(12,:,ii))); hold on; 
% end
% subplot(2,1,2); plot(db(sum(vibration_spectrum(12,:,:), 3)  )); 
% 
lowpass(sum(corrected_phases(8, :, :),3),  1500, 14925);
% spec_of_interest = sum(vibration_spectrum(8, :, :),3);
% noisesamples = abs(spec_of_interest(800:1024));
% estmean = mean(noisesamples);
% eststd = std(noisesamples);
% figure; plot([0:1023]/params.chirpCycleTime/params.opDopplerFFTSize,...
%     db(spec_of_interest(1:1024))); hold on; plot([1 14900/2], db(2*eststd+ [estmean estmean]));
% findpeaks(abs(spec_of_interest(1:400)), 'MinPeakHeight', 2*eststd+estmean)


% % filtering 
% for frame_idx = 9:14
%     lo = 1 + floor(70 / (1/params.chirpCycleTime/params.opDopplerFFTSize));
%     hi = 1 + ceil(160 / (1/params.chirpCycleTime/params.opDopplerFFTSize));
%     [M, I] = max(abs(vibration_spectrum(frame_idx, lo:hi)));
%     fundamental = (lo + (I-1) -1) * (1/params.chirpCycleTime/params.opDopplerFFTSize);
%     % init
%     frame_phases = corrected_phases(frame_idx, :); 
%     figure; lowpass(frame_phases, 2000, 14925);
% %     eq_phases = frame_phases; 
% %     eq_db = [0 20 27 32 34 37 40];
% %     for ii=1:5
% %         bpFilt2 = designfilt('bandpassfir','FilterOrder',2000, 'CutoffFrequency1',ii*fundamental-50,...
% %             'CutoffFrequency2',ii*fundamental+50, 'SampleRate',1/params.chirpCycleTime);
% %         phases2 = conv(bpFilt2.Coefficients, frame_phases); phases2 = phases2(1+1000:end-1000);
% %         eq_phases = eq_phases + 10^(eq_db(ii)/20) * phases2 - phases2;
% %     end
% %     [v, c, I] = phases_to_spectrum(eq_phases, params, polyFitOrder, 1);
% %     vibration_spectrum(frame_idx, :) = v;
% %     corrected_phases(frame_idx, :) =c;
% end

% spg_plot = vibration_spectrum(:, 1:params.opDopplerFFTSize/2).';
% figure; imagesc(abs(spg_plot)); colorbar;
% imagesc(db(spg_plot./max(max(abs(spg_plot)))));


for ii=1:4
    if ii==4
        y_spec = reshape(sum(corrected_phases, 3).', [], 1);
    else
         figure(221); subplot(2,2,ii); imagesc(abs(vibration_spectrum(:, 1:params.opDopplerFFTSize/2, ii).')); colorbar;
        y_spec = reshape(corrected_phases(:,:,ii).', [], 1);
    end
    Nx = length(y_spec);
    nsc = 1024; %params.numChirps; %floor(Nx/16);
    nov = floor(nsc/2);
    nff = max(8192,2^nextpow2(nsc));
    fs = 1/params.chirpCycleTime;
    figure(222);subplot(2,2,ii); spectrogram(y_spec,hamming(nsc),nov,nff, fs, 'yaxis');%,'MinThreshold',-98); 
    ylim([0 1.5]);
end
% figure;  imagesc(abs(s(:, 1:2:end)));

%%
for jj=1:1%params.numFrames
    phases1_fullband = zeros(1, params.numChirps);
    phases1_halfband = zeros(1, params.numChirps);
    phases2_fullband = zeros(1, params.numChirps);
    phases2_halfband = zeros(1, params.numChirps);
    for ii=1:params.numChirps
        adc_data = double(squeeze(datacube1.adcdata{jj}(ii,1,:))).';
        rangefft = fft(adc_data, params.opRangeFFTSize);
        phases1_fullband(ii) = angle(rangefft(peak));
        
        rangefft = fft(adc_data(1:1000), params.opRangeFFTSize);
        phases1_halfband(ii) = angle(rangefft(peak));
        
        adc_data = double(squeeze(datacube2.adcdata{jj}(ii,1,:))).';
        rangefft = fft(adc_data, params.opRangeFFTSize);
        phases2_fullband(ii) = angle(rangefft(peak));
        
        rangefft = fft(adc_data(10+[1:1000]), params.opRangeFFTSize);
        phases2_halfband(ii) = angle(rangefft(peak));
    end
    polyFitOrder = 9;
    [vibration_spectrum1_fullband, I] = phases_to_spectrum(phases1_fullband, params, polyFitOrder, 1);
    [vibration_spectrum1_halfband, I] = phases_to_spectrum(phases1_halfband(1:1000), params, polyFitOrder, 1);
    [vibration_spectrum2_fullband, I] = phases_to_spectrum(phases2_fullband, params, polyFitOrder, 1);
    [vibration_spectrum2_halfband, I] = phases_to_spectrum(phases2_halfband(14+[1:1000]), params, polyFitOrder, 1);
end
angle([vibration_spectrum1_fullband(I) vibration_spectrum1_halfband(I) ...
    vibration_spectrum2_fullband(I) vibration_spectrum2_halfband(I)])

[vibration_spectrum_sum, I] = phases_to_spectrum(...
    phases1_halfband(1:1000)+phases2_halfband(14+[1:1000]), params, polyFitOrder, 1);


% I = floor(1000/(1/params.chirpCycleTime/params.opDopplerFFTSize));
% abs(vibration_spectrum(I-1:I+1))
% noise_sample_idx = floor(params.opDopplerFFTSize/2*0.8):params.opDopplerFFTSize/2;
% db(abs(vibration_spectrum(I))/mean(abs(vibration_spectrum(noise_sample_idx))))


% % doppler FFT
% if plot_doppler_fft
%     rangeFFT_size = length(radarCube.data{selected_frame}(selected_chirp, selected_rxchain,:));%2^nextpow2(rfParams.numRangeBins);
%     dopplerFFT_size = 2^nextpow2(rfParams.numDopplerBins);
%     rangefft_fmcw_frame = radarCube.data{1, selected_frame};
%     rangefft_chirps = squeeze(rangefft_fmcw_frame(:,selected_rxchain,:));
%     dopplerfft = fft(rangefft_chirps, dopplerFFT_size, 1); % fft across chirps
%     cf = cf + 1; figure(cf);
%     imagesc(db(fftshift(dopplerfft, 1)));
%     xlabel('range (m)');
%     ylabel('velocity (m/s)');
%     range_bins = [0:rangeFFT_size/8:dims.numRangeBins-1];
%     ranges = rfParams.rangeResolutionsInMeters*range_bins;
%     xticks(range_bins);
%     xticklabels(arrayfun(@num2str, round(ranges,2), 'UniformOutput', false));
%     velocity_bins = [-dopplerFFT_size/2:dopplerFFT_size/8:dopplerFFT_size/2-1];
%     velocities = rfParams.dopplerResolutionMps*velocity_bins;
%     yticks(velocity_bins+dopplerFFT_size/2+1);
%     yticklabels(arrayfun(@num2str, round(velocities,2), 'UniformOutput', false));
%     colorbar;
% end

% phase across chirps
% if plot_phase_across_chirps
%     I=20; %object_dist = rfParams.rangeResolutionsInMeters * (I-1); %0.93m
% %     I=37; object_dist = rfParams.rangeResolutionsInMeters * (I-1); %1.547m
%     cf = cf + 1; figure(cf);
%     for x = 1:min(6, dims.numFrames)
%         subplot(2,3,x);
%         adc_data = double(squeeze(datacube.adcdata{x}(ii,1,:))).';
%         rangefft = fft(adc_data, params.opRangeFFTSize);
%         phases(:, ii) = angle(rangefft(peak));
%         
%         
%         rangefft_fmcw_frame = radarCube.data{x};
%         plot(angle(rangefft_fmcw_frame(1:end,selected_rxchain,I))/pi*180);
%         xlabel('chirp number'); ylabel('phase (deg)');
%         %title(sprintf('frame%d dist=%.3f m', x, object_dist));
%     end
% end

%% vibration frequency estimation, FFT on phases across chirps
close all;
tic
freq_interval = params.sampleRate*1e6/params.opRangeFFTSize;
dist_interval = freq_interval*3e8/(2*params.freqSlope*1e12);
peak_index = 170; %574 speaker, 335 trashbag
fprintf('object distance %.4f m\n', (peak_index-1)*dist_interval);

selected_frame = 1;
selected_rxchain = 1;

frame = datacube.adcdata{selected_frame};
clear chirp_phases    
N=1024; 
numSubChirp = 2;
subChirpDiricWindowSize = floor(params.numSamplePerChirp/numSubChirp);

% params.subChirpSincWindowSize = subChirpDiricWindowSize;   
% params.subChirpSincs = exp(-1j*(params.w)*(params.subChirpSincWindowSize-1)/2).*...
%                     diric(params.w, params.subChirpSincWindowSize);

numSincs = 10;
perChirp_w_idx = zeros(N, numSincs);
perChirp_rangefft_idx = zeros(N, numSincs);
perChirp_fft = zeros(N, params.opRangeFFTSize);
perChirp_sincs = zeros(N, numSincs, params.opRangeFFTSize);
perChirp_phis = zeros(N, numSincs);

perSubChirp_fft = zeros(N, numSubChirp, params.opRangeFFTSize);
perSubChirp_sincs = zeros(N, numSubChirp, numSincs, params.opRangeFFTSize);
perSubChirp_phis = zeros(N, numSubChirp, numSincs);

% gpd = mean(grpdelay(mybandpass, params.numSamplePerChirp, params.sampleRate*1e6));
for jj = 1:N
    adc_data = double(squeeze(frame(jj,selected_rxchain,:))).';
    
    % filtering
%     r1 = real(adc_data); i1 = imag(adc_data);
%     r2 = conv(mybandpass.', r1); i2 = conv(mybandpass.', i1);
%     r2 = r2(gpd+1:end-gpd); i2 = i2(gpd+1:end-gpd);
%     adc_data = r2 + 1j*i2;
    
    fftsize = params.opRangeFFTSize;
    w = 2*pi*[0:fftsize-1]/fftsize;
    diricWindowSize = length(adc_data);

%     if jj==1
%         [out2] = em_algo(adc_data, params, numSincs, ...
%             diricWindowSize,[], 0);
%     else
%         [out2] = em_algo(adc_data, params, numSincs, ...
%             diricWindowSize,out2, 0);
%     end
%     
%     perChirp_sincs(jj, :,:) = out2.sincs;
%     perChirp_phis(jj, :) = out2.phis;
%     perChirp_w_idx(jj,:) = out2.w_idx;
%     perChirp_rangefft_idx(jj,:) = out2.rangefft_idx;
    perChirp_fft(jj, :) = fft(adc_data, params.opRangeFFTSize);
%     fprintf('[%d] leftover norm %.1f %%\n',jj,...
%     100*norm(perChirp_fft(jj, :)-sum(out2.sincs,1))/norm(perChirp_fft(jj, :)));

    for ii=1:numSubChirp
%           subchirp_adcdata = adc_data((ii-1)*32+[1:512-32*3]);

        subchirp_adcdata = adc_data((ii-1)*subChirpDiricWindowSize+[1:subChirpDiricWindowSize]);

%         centers = params.w(out2.w_idx).';
%         phase_adjust = 4*pi*out2.dists/3e8*((ii-1)*...
%             subChirpDiricWindowSize*params.freqSlope/params.sampleRate*1e6);
% 
%         pred_sincs = (1/numSubChirp)*(out2.amps.*exp(1j*out2.phis).*...
%             exp(1j*phase_adjust)).'.*...
%             exp(-1j*(w-centers)*(subChirpDiricWindowSize-1)/2).*...
%             diric(w - centers, subChirpDiricWindowSize);
%     
%         out3 = out2;
%         out3.amps = (1/numSubChirp)*out2.amps;
%         out3.sincs = pred_sincs;
% 
%         [out3] = em_algo(subchirp_adcdata, params, numSincs,...
%             subChirpDiricWindowSize,out3,0);
%         
%         perSubChirp_sincs(jj,ii,:,:) = out3.sincs;
%         perSubChirp_phis(jj, ii,:) = out3.phis;
        perSubChirp_fft(jj,ii, :) = fft(subchirp_adcdata, fftsize);
    end
    
%     figure(100);
%     subplot(5,1,jj); plot(abs(perChirp_fft(jj, :))); hold on;
%     plot(abs(sum(sincs(:, :),1))); xlim([0 600]);
end
toc
%%
sum(abs(perChirp_w_idx(:,2)-285262))
perChirp_w_idx(1:10,1:5)
perChirp_rangefft_idx(1,1:end)
(perChirp_rangefft_idx(1,1:end)-1)*dist_interval

figure; plot(perChirp_w_idx(:,6));
% 1. chirp fft
chirp_phases = angle(perChirp_fft(:, [22]).'); 
chirp_phases = angle(perChirp_fft(:, [172]).'); 
chirp_phases = angle(perChirp_fft(:, [58]).'); 

% 2. chirp sinc
chirp_phases = perChirp_phis(:,:).'; 
chirpleftover = perChirp_fft - squeeze(sum(perChirp_sincs(:,:,:), 2));
chirp_phases = angle(perChirp_sincs(:,[1:10],286).'); 
chirp_phases = angle(squeeze(perChirp_sincs(:,[3],...
    perChirp_rangefft_idx(1,1:end))).'); 
chirp_phases = angle(chirpleftover(:,294).'); 

chirp_phases = angle(sum(perChirp_sincs(:,[2 6 10],286), 2).'); 

chirp_phases = angle((sum(perChirp_sincs(:,[2 6 10],286), 2)+...
    chirpleftover(:,286)).'); 
chirp_phases = angle(sum(perChirp_sincs(:,[3 8],172), 2).'); 

chirp_phases = angle(sum(perChirp_sincs(:,[3 10],172), 2).'); 

% chirp_phases = perChirp_phis(:,:,2).';
% 3. subchirp fft
chirp_phases = angle(perSubChirp_fft(:, :, 22).'); 

% 4. subchirp sinc
subChirpleftover = perSubChirp_fft - squeeze(sum(perSubChirp_sincs, 3)); 
chirp_phases = angle(perSubChirp_sincs(:,:,10,286).' + subChirpleftover(:,:,286).'); %17.9db
chirp_phases = angle(squeeze(perSubChirp_sincs(:,2,[1:15],286)).');
chirp_phases = squeeze(perSubChirp_phis(:,1,:)).';

subChirpleftover = perSubChirp_fft - squeeze(sum(perSubChirp_sincs(:,:,[12],:),3)); 
chirp_phases = angle(subChirpleftover(:, :, 286).'); 

chirp_phases = angle(sum(perSubChirp_sincs(:,:,[1],285),3).');
chirp_phases = angle((sum(perSubChirp_sincs(:,:,[3 6 10],286),3)+...
   subChirpleftover(:,:,286)).');
chirp_phases = angle(sum(perSubChirp_sincs(:,:,[3 8],172),3).');

% chirp_phases(2,:) = mod(chirp_phases(2,:) - 2, 2*pi);
%%
close all;
polyFitOrder = 9;
[vibration_spectrum, I] = phases_to_spectrum(chirp_phases, params, polyFitOrder, 1);
% onelongchirp_noise = db(mean(abs(vibration_spectrum(:,fftsize/2-500:fftsize/2)), 2));
% onelongchirp_snr= db(abs(vibration_spectrum(:,I))) - onelongchirp_noise;
% fprintf('snr %.1f dB\n', onelongchirp_snr);
% % 
abs(vibration_spectrum(:, 31))
angle(vibration_spectrum(:, 31))

sum_spec = sum(vibration_spectrum, 1);
figure; plot(abs(sum_spec));
noise_sample_idx = floor(params.opDopplerFFTSize/2*0.8):params.opDopplerFFTSize/2;
db(abs(sum_spec(31))/mean(abs(sum_spec(noise_sample_idx))))
db(abs(sum_spec(42))/mean(abs(sum_spec(1024-500:1024))))

abs(sum_spec(826))
figure; plot(1:params.opDopplerFFTSize, ...
    angle(vibration_spectrum(1,:)./vibration_spectrum(2:end,:)));
legend('subchirp2', 'subchirp3', 'subchirp4'); 
ylabel('phase difference'); xlabel('vibration bin');

figure;
for ii=1:4
    plot(perSubChirp_fft(:, ii, 405), 'o'); hold on;
end
% for ii=1:size(chirp_phases, 1)
%     t = 1:length(chirp_phases(ii,:));
%     [p, s, mu] = polyfit(t, chirp_phases(ii,:), polyFitOrder);
%     trend(ii,:) = polyval(p, t, [], mu);
%     corrected_phases(ii,:) = chirp_phases(ii,:) - trend(ii,:);
% %     fftout(ii,:) = fft(corrected_phases(ii,:), params.opDopplerFFTSize);
% end
% figure;
% dfft = fft(corrected_phases, params.opDopplerFFTSize, 2);
% imagesc(abs(dfft)./max(abs(dfft),[], 2));
% colorbar;
% figure; plot(abs(dfft(286,:)));
%%
% 
% radarCube_ChirpCorrected = chirp_phase_correction(radarCube, peak_index, 0);
% phases_original = angle(radarCube.data{1}(:,1,peak_index));
% phases_corrected = angle(radarCube_ChirpCorrected.data{1}(:,1,peak_index));

numSubChirps = 2;
subChirpStrideStep = floor(params.numSamplePerChirp/numSubChirps);
subChirpSamples = floor(params.numSamplePerChirp/numSubChirps);
clear chirp_subchirp_phases chirp_phases corrected_phases spectrum subchirp_snr subchirp_noise
for jj = 1:params.numChirps
    adc_data = squeeze(frame(jj,selected_rxchain,:)).';
    for hh = 1:numSubChirps
        start_index = (hh-1)*subChirpStrideStep+1;
        fft_output = fft(adc_data(start_index:start_index+subChirpSamples-1), params.opRangeFFTSize);     
        chirp_subchirp_phases(jj, hh) = angle(fft_output(peak_index));
        
        if jj==1
            figure(5);
            subplot(numSubChirps+1, 1, 1);
            plot(dist_interval*[0:params.opRangeFFTSize-1], abs(fft(adc_data, params.opRangeFFTSize)));
            title('Original long chirp Range FFT'); xlim([0 5]);
            subplot(numSubChirps+1, 1, hh+1);
            plot(dist_interval*[0:params.opRangeFFTSize-1], abs(fft_output));
            title(sprintf('Subchirp %d Range FFT', hh)); xlim([0 5]);
        end
    end
end

for ii=1:numSubChirps
    chirp_phases = chirp_subchirp_phases(:, ii).';
    t = 1:length(chirp_phases);
    [p, s, mu] = polyfit(t, chirp_phases, 9);
    trend = polyval(p, t, [], mu);
    corrected_phases(ii,:) = chirp_phases - trend;
    spectrum(ii, :) = fft(corrected_phases(ii,:), params.opDopplerFFTSize);
    
    figure(7);
    subplot(numSubChirps, 2, (ii-1)*2+1); plot(chirp_phases); hold on; plot(trend);
    subplot(numSubChirps, 2, (ii-1)*2+2); plot(corrected_phases(ii,:));
%     highpass(corrected_phases(4,:), [100], 1/params.chirpCycleTime)
    
    figure(6);
    subplot(numSubChirps, 1, ii);
    Fs = 1/params.chirpCycleTime;
    plot(Fs/1e3/params.opDopplerFFTSize.*[0:params.opDopplerFFTSize-1], abs(spectrum(ii, :)));
    xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');
end
[M, I] = max(abs(spectrum(:,1:params.opDopplerFFTSize/2)), [], 2);
mode_I = mode(I);
fprintf('Subchirps majority frequency peak at %.3f kHz Index=%d', Fs/1e3/params.opDopplerFFTSize*(mode_I-1), mode_I);
disp(I.');
subchirp_noise = db(mean(abs(spectrum(:, params.opDopplerFFTSize/2-500:params.opDopplerFFTSize/2)), 2)).';
subchirp_snr = db(abs(spectrum(:, mode_I))).' - subchirp_noise;


fprintf('subchirp SNR (dB) : '); disp(subchirp_snr);
fprintf('subchirp noise (dB) : '); disp(subchirp_noise);
fprintf('onelongchirp SNR (dB) : '); disp(onelongchirp_snr);
fprintf('onelongchirp noise (dB) : '); disp(onelongchirp_noise);

phases_combined = sum(corrected_phases(1:end,:), 1);
spectrum_combined = fft(phases_combined, params.opDopplerFFTSize);
[M, I] = max(abs(spectrum_combined(1:params.opDopplerFFTSize/2)));
combined_noise = db(mean(abs(spectrum_combined(params.opDopplerFFTSize/2-500:params.opDopplerFFTSize/2))));
combined_snr = db(abs(spectrum_combined(I))) - combined_noise;
fprintf('combined subchirps SNR (dB) : '); disp(combined_snr);
fprintf('combined subchirps (dB) : '); disp(combined_noise);
fprintf('SNR improvement? (dB) %.3f --> %.3f\n', onelongchirp_snr,combined_snr);

figure;
plot(abs(a1)); hold on; plot(abs(spectrum_combined));


candidate = [mode_I 200 300 400 500 600];
figure(50);
subplot(1,2,1);
for ss=candidate
    plot(1:numSubChirps+1, angle([spectrum(:, ss); spectrum_combined(ss);]), '^-'); hold on;
    plot(numSubChirps+1, angle(spectrum_combined(ss)), 'k^'); 
    xlabel('subchirps'); title('phase');
end
% plot(bad_subchirps, angle(spectrum(bad_subchirps, mode(I))), 'ro');
subplot(1,2,2);
for ss=candidate
    plot(1:numSubChirps+1, abs([spectrum(:, ss); spectrum_combined(ss);]), '^-'); hold on;
    plot(numSubChirps+1, abs(spectrum_combined(ss)), 'k^'); 
    xlabel('subchirps'); title('amplitude');
end



%% experimental scripts
% combine multiple subchirps
Fs = 1/45e-6;
clear radarCube_SubChirp
radarCube_SubChirp = get_subchirps(radarCube, 128, 128);
frame1 = radarCube_SubChirp.data{1};
% verify all subchirps have the same peak
for ii=1:radarCube_SubChirp.dim.numSubChirps
    [M, I] = max(abs(squeeze(frame1(1,ii,1,:))));
    maxPos(ii) = I;
end
figure;
subplot(2,1,1); plot(maxPos);
subplot(2,1,2); plot(abs(squeeze(frame1(1,1,1,:))));

peak_index = 10;
fprintf('object distance %.4f m\n', (peak_index-1)*radarCube_SubChirp.rfParams.rangeResolutionsInMeters);

chirp_subchirps = radarCube_SubChirp.data{1}(:,:,1,peak_index);
clear raw_phases3 corrected_phases3
for ss=1:radarCube_SubChirp.dim.numSubChirps
    chirp_phases =  mod(angle(chirp_subchirps(:,ss)), 2*pi).';
    raw_phases3(ss, :) = chirp_phases;
    t = 1:length(chirp_phases);
    [p, s, mu] = polyfit(t, chirp_phases, 9);
    trend = polyval(p, t, [], mu);
    corrected_phases3(ss, :) = chirp_phases - trend;
end
clear spectrum spectrum_normalized spectrum_combined spectrum_combined_normalized
spectrum = fft(corrected_phases3, 2048, 2);
spectrum_normalized = fft(corrected_phases3./max(corrected_phases3, [], 2), 2048, 2);
[M, I] = max(abs(spectrum), [], 2);
mode(I)

% plot bad (noisy) subchirps
bad_subchirps = find(abs(I-mode(I))>1);
L = length(bad_subchirps);
figure;
for ss=1:L
    subplot(ceil(L/5), 5, ss); plot(abs(spectrum(bad_subchirps(ss), :)));
    title(sprintf('subchirp %d', bad_subchirps(ss)));
end

% plot phases/amplitudes of subchirps
candidate = [mode(I) 82:300];
spectrum_combined = fft(sum(corrected_phases3, 1));
corrected_phases3_normalized = corrected_phases3./max(corrected_phases3, [], 2);
spectrum_combined_normalized = fft(sum(corrected_phases3_normalized, 1));
N_subchirps = radarCube_SubChirp.dim.numSubChirps;
figure(50);
subplot(1,2,1);
for ss=candidate
    plot(1:N_subchirps+2, angle([spectrum(:, ss); spectrum_combined(ss); spectrum_combined_normalized(ss)]), '^-'); hold on;
    plot(N_subchirps+1, angle(spectrum_combined(ss)), 'k^'); hold on;
    plot(N_subchirps+2, angle(spectrum_combined_normalized(ss)), 'k^'); 
    xlabel('subchirps'); title('phase');
end
% plot(bad_subchirps, angle(spectrum(bad_subchirps, mode(I))), 'ro');
subplot(1,2,2);
for ss=candidate
    plot(1:N_subchirps+2, abs([spectrum(:, ss); spectrum_combined(ss); spectrum_combined_normalized(ss)]), '^-'); hold on;
    plot(N_subchirps+1, abs(spectrum_combined(ss)), 'k^'); hold on;
    plot(N_subchirps+2, abs(spectrum_combined_normalized(ss)), 'k^'); 
    xlabel('subchirps'); title('amplitude');
end
% plot(bad_subchirps, abs(spectrum(bad_subchirps, mode(I))), 'ro');
figure;
chirp1024_fft = abs(a1);
plot(db(chirp1024_fft(candidate(1))./chirp1024_fft(candidate(2:end)) )); hold on;
plot(db(abs(spectrum_combined(candidate(1)))./abs(spectrum_combined(candidate(2:end))) )); hold on;
plot(db(abs(spectrum_combined_normalized(candidate(1)))./abs(spectrum_combined_normalized(candidate(2:end)) ))); hold on;

x1 = db(chirp1024_fft(candidate(1))./mean(chirp1024_fft(100:300)) );
x2 = db(abs(spectrum_combined(candidate(1)))./mean(abs(spectrum_combined(100:300))) );
x3 = db(abs(spectrum_combined_normalized(candidate(1)))./mean(abs(spectrum_combined_normalized(100:300))) );
[x1 x2 x3]
%%


angle([spectrum(:,250) spectrum_normalized(:, 250)])
angle([spectrum_combined(250) spectrum_combined_normalized(250)])
abs([spectrum(:,250); spectrum_combined(250); spectrum_combined_normalized(250)])

figure(60);
for ss=1:4
    phases = corrected_phases3(ss, :);
    subplot(4,2,(ss-1)*2+1);
    plot(phases);
    subplot(4,2,(ss-1)*2+2);
    plot(Fs/1e3/2048.*[0:2047], abs(fft(phases)));
    xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');
end
combine = sum(corrected_phases3, 1);
figure(61); plot(Fs/1e3/2048.*[0:2047], abs(fft(combine)));
xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');

figure(62);
for ss=1:4
    phases = corrected_phases3(ss, :);
    subplot(4,2,(ss-1)*2+1);
    plot(phases./max(phases));
    subplot(4,2,(ss-1)*2+2);
    plot(Fs/1e3/2048.*[0:2047], abs(fft(phases./max(phases))));
    xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');
end
combine = sum(corrected_phases3./max(corrected_phases3, [], 2), 1);
figure(63); plot(Fs/1e3/2048.*[0:2047], abs(fft(combine)));
xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');


figure;
for ss=1:4
%     subplot(4,1,ss); plot(raw_phases3(ss, : ));
%     plot(raw_phases3(ss, : )); hold on;
    plot(corrected_phases3(ss, 1:256 )); hold on;
end
raw_phases3(:, 100:103)

% see how subchirps combine
candidate = [1:4];
L = length(candidate);
figure;
for ss=1:L
    subplot(ceil(L/5), 5, ss); plot(abs(spectrum(candidate(ss), :)));
    title(sprintf('subchirp %d', candidate(ss)));
end
combine = sum(corrected_phases3(candidate, :), 1);
subplot(ceil(L/5), 5, 5); plot(abs(fft(combine)));title('combined');


% radarCube_SubChirpCorrected = subchirp_phase_correction(radarCube_SubChirp, peak_index, 'characterization');
% radarCube_SubChirpCorrected = subchirp_phase_correction(radarCube_SubChirp, peak_index, 'fitted');
% radarCube_ChirpCorrected = chirp_phase_correction(radarCube_SubChirpCorrected, peak_index, 1);

%%
% phases across chirps by rx chains
selected_frame = 3;
Fs = 1/45e-6;
rangefft_fmcw_frame = radarCube.data{selected_frame};
dopplerFFT_size = 2^nextpow2(rfParams.numDopplerBins);
for rxchain=1:4
    chirps = rangefft_fmcw_frame(1:end,rxchain,I).';
    y = angle(chirps);
    t = 1:length(y);
    [p, s, mu] = polyfit(t, y, 9);
    y_trend = polyval(p, t, [], mu);
    dt_y(rxchain,:) = y - y_trend;
    figure(4);
    subplot(3,4,rxchain); plot(abs(chirps)); 
    title(sprintf('rx%d amplitude', rxchain)); xlabel('chirps');
    subplot(3,4,rxchain+4); plot(y); hold on; plot(y_trend); 
    legend('phase samples', 'trend', 'Location', 'southeast'); title(sprintf('rx%d original phase', rxchain));
    subplot(3,4,rxchain+8); plot(dt_y(rxchain,:)); 
    title(sprintf('rx%d detrended phase', rxchain));
end
% doppler FFT, detect vibration frequency, across rx chains
for rxchain=1:4
    figure(5);
    % filtered_dt_y = lowpass(dt_y, 3000, Fs);
    subplot(2,4,rxchain); 
    plot(Fs/1e3/dopplerFFT_size.*[0:dopplerFFT_size-1], abs(fft(dt_y(rxchain,:), dopplerFFT_size)));
    xlim([0 Fs/1e3]); xlabel('Frequency (kHz)'); title(sprintf('rx%d %d-FFT', rxchain, dopplerFFT_size));
    subplot(2,4,rxchain+4); 
    plot(Fs/1e3/1024.*[0:1024-1], abs(fft(dt_y(rxchain, 1025:end), 1024)));
    xlim([0 Fs/1e3]); xlabel('Frequency (kHz)'); title(sprintf('rx%d %d-FFT', rxchain, 1024));
    % plot(abs(fft(filtered_dt_y, dopplerFFT_size)));
    % subplot(2,1,2); plot(abs(fft(chirps, dopplerFFT_size))); % this doesn't work
end
% doppler FFT, detect vibration frequency, across frames
figure(6);
idx = [17 22 37]; % door, trachbag, speaker
frames = min(6, dims.numFrames);
for frame=1:frames
    rangefft_fmcw_frame = radarCube.data{frame};
    for jj = 1:3
        subplot(frames,3,(frame-1)*3 + jj); 
        I = idx(jj);
        object_dist = rfParams.rangeResolutionsInMeters * (I-1);
        y = angle(rangefft_fmcw_frame(1:end,selected_rxchain,I).');
        t = 1:length(y);
        [p, s, mu] = polyfit(t, y, 9);
        y_trend = polyval(p, t, [], mu);
        dt_y = y - y_trend;
        plot(Fs/1e3/dopplerFFT_size.*[0:dopplerFFT_size-1], abs(fft(dt_y, dopplerFFT_size)));
        xlim([0 Fs/1e3]); xlabel('Frequency (kHz)'); title(sprintf('frame%d d=%.3fm', frame, object_dist));
    end
end

for ii=1:4
    rangefft_fmcw_frame = radarCube.data{ii};
    chirps = rangefft_fmcw_frame(1:end,1,37).';
    plot(chirps, 'o'); hold on;
end
% figure(7);
% spectrogram = zeros(256, 700);
% for ii=1:700
%     rangefft_fmcw_frame = radarCube.data{ii};
%     x = angle(rangefft_fmcw_frame(1:end,selected_rxchain,I));
%     x = detrend(x);
%     spectrogram(:, ii) = db(fft(x, 256));
% end
% imagesc(spectrogram);
% colorbar;
%% post processing -- multi-frame, single chirp per frame
selected_rxchain = 1;
rangefft_fmcw_frames = zeros(dims.numFrames, dims.numRangeBins);
for idx = 1 : dims.numFrames
    rangefft_fmcw_frames(idx, :) = squeeze(radarCube.data{idx}(1,selected_rxchain,:));
end
[M, I] = max(abs(rangefft_fmcw_frames(1,:)));
figure(1);
subplot(2,2,1); plot(abs(rangefft_fmcw_frames(:,I)));
subplot(2,2,2); plot(angle(rangefft_fmcw_frames(:,I))); xlabel('chirps'); ylabel('phase (rad)');
subplot(2,2,3); plot(rfParams.rangeResolutionsInMeters*[0:dims.numRangeBins-1], db(rangefft_fmcw_frames(1,:))); xlabel('range (m)');
dopplerfft = fft(rangefft_fmcw_frames(:,I));
subplot(2,2,4); plot(abs(dopplerfft));
