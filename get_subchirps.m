function [newradarCube] = get_subchirps(radarCube, new_rangefft_size, subChirpStrideStep)
    N_frames = 1; %length(radarCube.data);
    N_chirps = radarCube.dim.numChirps;
    N_rxchains = 1; %radarCube.dim.numRxChan;

    newradarCube.dim = radarCube.dim;
    newradarCube.rfParams = radarCube.rfParams;
    numSubChirps = ceil( (radarCube.dim.numRangeBins-new_rangefft_size+1)/subChirpStrideStep);
    
    newradarCube.dim.numSubChirps = numSubChirps;
    newradarCube.dim.numRangeBins = new_rangefft_size;
    newradarCube.dim.subChirpStrideStep = subChirpStrideStep;
    
    newradarCube.rfParams.numRangeBins = new_rangefft_size;
    newradarCube.rfParams.bandwidth = radarCube.rfParams.freqSlope * (1/radarCube.rfParams.sampleRate)*1e6*new_rangefft_size/1e9;
    newradarCube.rfParams.rangeResolutionsInMeters = 3e8 /2/(newradarCube.rfParams.bandwidth*1e9);
    
    fprintf('Range FFT size: %d --> %d\n', radarCube.dim.numRangeBins, new_rangefft_size);
    fprintf('Range resolution(m): %.3f --> %.3f\n', radarCube.rfParams.rangeResolutionsInMeters,...
        newradarCube.rfParams.rangeResolutionsInMeters);
    
    for ii = 1:N_frames
        frame = radarCube.data{ii};
        newframe = zeros(radarCube.dim.numChirps,numSubChirps,radarCube.dim.numRxChan,new_rangefft_size);
        for jj = 1:N_chirps
            for kk = 1:N_rxchains
                chirp_at_rx = squeeze(frame(jj,kk,:));
                raw_adc_data = ifft(chirp_at_rx);
                for hh = 1:numSubChirps
                    start_index = (hh-1)*subChirpStrideStep+1;
                    newframe(jj, hh, kk, :) = fft(raw_adc_data(start_index:start_index+new_rangefft_size-1));
                end
            end
        end
        newradarCube.data{ii} = newframe;
    end
end