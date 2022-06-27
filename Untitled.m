tic
clear all; close all;
savefile = 0;
loadfromfile = 0;
foldername = 'experiment_data/layschip/radarcube';
filename = '0_radar0.5m_source1.0m_profile3_voice0_2';
wavpath =  'experiment_data/layschip/wav_noeq';
offset = 6;
params = read_from_json('profile3.mmwave.json');
numFrames = 7;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.numFrames = numFrames;
if loadfromfile
    load(sprintf('%s/%s.mat', foldername, filename));
else
    numBinFiles = ceil((4*params.numSamplePerChirp*4)*params.numChirps*numFrames/1024^3);
    adc_data =  read_from_binfile(numBinFiles, numFrames, params.numChirps, params.numSamplePerChirp);
    % adc_data (chirpsamples, chirps, frames)
    datacube.adcdata = adc_data;
    datacube.params = params;
end

% range fft
[fftout, I] = rangeFFT(datacube.adcdata(:,1,1), datacube.params);
% extract phases across chirps
rangecube = fft(datacube.adcdata, datacube.params.opRangeFFTSize, 1);
y = angle(reshape(rangecube(I,:,:), [], 1)); % TODO: piecewise phase correction
y2 = zeros(size(y));
segment_size = 1000;
% phase correction
for ii=1:datacube.params.numFrames
    for jj=1:segment_size:datacube.params.numChirps
        start_idx = (ii-1)*datacube.params.numChirps + jj;
        end_idx = min(start_idx+segment_size-1, ii*datacube.params.numChirps); % don't cross chirp border
        t = [1:end_idx-start_idx+1].';
        [p, s, mu] = polyfit(t, y(start_idx:end_idx), 9);
        trend= polyval(p, t, [], mu);
        y2(start_idx:end_idx) = y(start_idx:end_idx) - trend;
    end
end
% figure; plot(y);  figure; plot(y2);   sound(0.01*y2/max(y2), 14925)

% spectrogram
nsc = 1024; %params.numChirps; %floor(Nx/16);
nov = floor(nsc/2);
nff = max(4096,2^nextpow2(nsc));
fs = 1/datacube.params.chirpCycleTime;
figure(100); spectrogram(y2,hamming(nsc),nov,nff, fs, 'yaxis','MinThreshold',-150); 
ylim([0 2]);

if savefile
    if ~exist(foldername, 'dir')
        mkdir(foldername);
    end
    save(sprintf('%s/%s.mat', foldername, filename), 'datacube', '-v7.3');
    savefig(gcf, sprintf('%s/%s.fig', foldername, filename));
end

toc
% sound(y2/max(y2), 1/datacube.params.chirpCycleTime);
% clear sound


%%%% cut out each utterance
clear startidx endidx
yyy = fft(reshape(y2(1:end-mod(length(y2),nsc)), nsc, []),nff);
% figure(101); imagesc(db(yyy)); 
zzz = max(abs(yyy), [],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = min(zzz) + 0.05*(max(zzz) - min(zzz));
utterance_period_threshold = 4; % in units of nsc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure; plot(zzz); hold on; plot([1 length(zzz)], threshold*ones(1,2));
truthvalue = zzz > threshold;
edge = truthvalue(1:end-1) - truthvalue(2:end);
posedge = 0;
cnt=0;
jj=1;
for ii=1:size(zzz, 2)-1
    if edge(ii) ==-1
        posedge = 1;
        cnt=0;
        startidx(jj) = max(ii-2, 1);
    elseif edge(ii)==1
        posedge = 0;
        if  cnt >utterance_period_threshold
            endidx(jj) = min(ii+2, size(zzz, 2));
            jj = jj +1;
        end
        cnt = 0;
    else 
        if posedge
            cnt = cnt + 1;
        end
    end
end

% rows = (max(endidx - startidx(1:length(endidx))) +1 ) * nsc;
% extracted_phases = zeros(rows, length(endidx));

if ~exist(wavpath, 'dir')
    mkdir(wavpath);
end

figure(102); subplot(2, length(endidx), length(endidx)+[1:length(endidx)]); plot(zzz); hold on ;
for ii=1:length(endidx)
    figure(102); subplot(2, length(endidx), length(endidx)+[1:length(endidx)]); 
    plot([startidx(ii) endidx(ii)], 0.15*max(zzz)*ones(1,2) ); hold on;
    figure(102); subplot(2,length(endidx),ii);
    phases = y2((startidx(ii)-1)*nsc+1:endidx(ii)*nsc);
    spectrogram(phases,hamming(nsc),nov,nff, fs, 'yaxis'); ylim([0 2]); 
%     save(sprintf('%s/%s%d.mat', wavpath, filename(1:end-1), offset+ii), 'phases');
%     audiowrite(sprintf('%s/%s%d.wav',wavpath, filename(1:end-1), ii+offset), phases, floor(1/datacube.params.chirpCycleTime));
end