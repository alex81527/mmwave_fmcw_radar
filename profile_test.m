close all;
tic
freq_interval = params.sampleRate*1e6/params.opRangeFFTSize;
dist_interval = freq_interval*3e8/(2*params.freqSlope*1e12);
peak_index = 172; %574 speaker, 335 trashbag
fprintf('object distance %.4f m\n', (peak_index-1)*dist_interval);

selected_frame = 1;
selected_rxchain = 1;

frame = datacube.adcdata{selected_frame};
clear chirp_phases    
N=5; 
numSubChirp = 2;
subChirpDiricWindowSize = floor(params.numSamplePerChirp/numSubChirp);

params.subChirpSincWindowSize = subChirpDiricWindowSize;   
params.subChirpSincs = exp(-1j*(params.w)*(params.subChirpSincWindowSize-1)/2).*...
                    diric(params.w, params.subChirpSincWindowSize);

numSincs = 15;
perChirp_w_idx = zeros(N, numSincs);
perChirp_fft = zeros(N, params.opRangeFFTSize);
perSubChirp_fft = zeros(N, numSubChirp, params.opRangeFFTSize);
perChirp_sincs = zeros(N, numSincs, params.opRangeFFTSize);
perSubChirp_sincs = zeros(N, numSubChirp, numSincs, params.opRangeFFTSize);
perSubChirp_phis = zeros(N, numSubChirp, numSincs);
for jj = 1:N
    adc_data = double(squeeze(frame(jj,selected_rxchain,:))).';
    fftsize = params.opRangeFFTSize;
    w = 2*pi*[0:fftsize-1]/fftsize;
    diricWindowSize = length(adc_data);

    if 1==1
        [out2] = em_algo(adc_data, params, numSincs, ...
            diricWindowSize,[], 0);
%     else
%         [out2] = em_algo(adc_data, params, numSincs, ...
%             diricWindowSize,out2, 0);
    end
    
    perChirp_sincs(jj, :,:) = out2.sincs;
    perChirp_w_idx(jj,:) = out2.w_idx;
    perChirp_fft(jj, :) = fft(adc_data, params.opRangeFFTSize);
    

    for ii=1:numSubChirp
        subchirp_adcdata = adc_data((ii-1)*subChirpDiricWindowSize+[1:subChirpDiricWindowSize]);

        centers = params.w(out2.w_idx).';
        phase_adjust = 4*pi*out2.dists/3e8*((ii-1)*...
            subChirpDiricWindowSize*params.freqSlope/params.sampleRate*1e6);

        pred_sincs = (1/numSubChirp)*(out2.amps.*exp(1j*out2.phis).*...
            exp(1j*phase_adjust)).'.*...
            exp(-1j*(w-centers)*(subChirpDiricWindowSize-1)/2).*...
            diric(w - centers, subChirpDiricWindowSize);
    
        out3 = out2;
        out3.amps = (1/numSubChirp)*out2.amps;
        out3.sincs = pred_sincs;

        [out3] = em_algo(subchirp_adcdata, params, numSincs,...
            subChirpDiricWindowSize,out3,0);
        
        perSubChirp_sincs(jj,ii,:,:) = out3.sincs;
        perSubChirp_phis(jj, ii,:) = out3.phis;
        perSubChirp_fft(jj,ii, :) = fft(subchirp_adcdata, fftsize);
    end
%     figure(100);
%     subplot(5,1,jj); plot(abs(perChirp_fft(jj, :))); hold on;
%     plot(abs(sum(sincs(:, :),1))); xlim([0 600]);
end
toc% 0.26/4
