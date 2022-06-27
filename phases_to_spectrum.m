function [fftout, corrected_phases, I] = phases_to_spectrum(chirp_phases, params, polyFitOrder, plotting)
    nChirps = size(chirp_phases, 1);
    chirpLen = size(chirp_phases, 2);
    I = zeros(nChirps, 1);
    
    chirp_phases = unwrap(angle(exp(1j*chirp_phases)),[], 2);
    for ii=1:nChirps
        t = 1:length(chirp_phases(ii,:));
        [p, s, mu] = polyfit(t, chirp_phases(ii,:), polyFitOrder);
        trend(ii,:) = polyval(p, t, [], mu);
        corrected_phases(ii,:) = chirp_phases(ii,:) - trend(ii,:);
        fftout(ii,:) = fft(corrected_phases(ii,:), params.opDopplerFFTSize);
    end
    
    if plotting 
        figure(101);
        for ii=1:nChirps
            subplot(nChirps,2,(ii-1)*2+1); 
            plot(chirp_phases(ii,:)); hold on; plot(trend(ii,:));
            subplot(nChirps,2,(ii-1)*2+2); 
            plot(corrected_phases(ii,:)); title('corrected phases');
        end
        figure(102);
        for ii=1:nChirps
            subplot(nChirps,1,ii); 
            Fs = 1/params.chirpCycleTime;
            plot(Fs/1e3/params.opDopplerFFTSize.*[0:params.opDopplerFFTSize-1], abs(fftout(ii,:)));
            xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');
            [M, I] = max(abs(fftout(ii, 1:params.opDopplerFFTSize/2)));
            
            noise_sample_idx = floor(params.opDopplerFFTSize/2*0.8):params.opDopplerFFTSize/2;
            noise = db(mean(abs(fftout(ii,noise_sample_idx))));
            snr = db(abs(fftout(ii,I))) - noise;
            fprintf('[%d] Frequency peak at %.3f kHz Index=%d snr=%.1f dB\n', ...
                ii,Fs/1e3/params.opDopplerFFTSize*(I-1), I, snr);
        end
    end
end