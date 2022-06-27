function [newradarCube] = chirp_phase_correction(radarCube, rangeBinIndex, correction_mode)
    N_frames = 1; %length(radarCube.data);
    N_chirps = radarCube.dim.numChirps;
    N_rxchains = 1; %radarCube.dim.numRxChan;
    
    if correction_mode == 1 % has subchirps
        for ii = 1:N_frames
            frame = radarCube.data{ii};
            for kk = 1:N_rxchains
                chirp_subchirps = squeeze(frame(:,:,kk,rangeBinIndex));%shape(N_chirps,N_subchirps)
%                 figure;
%                 for ss=1:4
%                     subplot(4,1,ss); plot(angle(chirp_subchirps(:,ss)));
%                 end
              %% 
                % 1
                chirps = radarCube.data{1}(:,1,37);
                chirp_phases =  mod(angle(chirps), 2*pi).';
                t = 1:length(chirp_phases);
                [p, s, mu] = polyfit(t, chirp_phases, 9);
                trend = polyval(p, t, [], mu);
                corrected_phases1 = chirp_phases - trend;
%                 frame(:,kk,rangeBinIndex) = abs(chirps_at_rx_rangebin).*exp(1j*corrected_phases);
                % 2
                chirp_subchirps = radarCube_SubChirpCorrected.data{1}(:,:,1,19);
                combined_chirps = sum(chirp_subchirps, 2);
                chirp_phases =  mod(angle(combined_chirps), 2*pi).';
                t = 1:length(chirp_phases);
                [p, s, mu] = polyfit(t, chirp_phases, 9);
                trend = polyval(p, t, [], mu);
                corrected_phases2 = chirp_phases - trend;
                % 3
                chirp_subchirps = radarCube_SubChirp.data{1}(:,:,1,4);
                for ss=1:radarCube_SubChirp.dim.numSubChirps
                    chirp_phases =  mod(angle(chirp_subchirps(:,ss)), 2*pi).';
                    t = 1:length(chirp_phases);
                    [p, s, mu] = polyfit(t, chirp_phases, 9);
                    trend = polyval(p, t, [], mu);
                    corrected_phases3(ss, :) = chirp_phases - trend;
                end
    
                spectrum = fft(corrected_phases3, 2048, 2);
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
                figure(50);
                plot(angle(spectrum(:, mode(I)))); hold on;
                plot(bad_subchirps, angle(spectrum(bad_subchirps, mode(I))), 'ro'); 
                hold on;
                plot(angle(spectrum(:, 150))); hold on ;
                plot(angle(spectrum(:, 180))); hold on ;
                plot(angle(spectrum(:, 280))); 
                
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
                
%                 y = bandpass(corrected_phases, [100 2000], Fs);    
                figure; 
                fftsize = 2048;
                plot(Fs/1e3/fftsize.*[0:fftsize-1], abs(fft(corrected_phases1))); 
                hold on; plot(Fs/1e3/fftsize.*[0:fftsize-1], abs(fft(combine3)));
                xlim([0 Fs/1e3]); xlabel('Frequency (kHz)');
              %%  
                figure;
                subplot(2,1,1); plot(corrected_phases1); %hold on; plot(trend);
                subplot(2,1,2); plot(corrected_phases2);     
                
            end
            radarCube.data{ii} = frame;
        end
    else 
        for ii = 1:N_frames
            frame = radarCube.data{ii};
            for kk = 1:N_rxchains
                chirps_at_rx_rangebin = squeeze(frame(:,kk,rangeBinIndex)).';
                chirp_phases =  angle(chirps_at_rx_rangebin);
                t = 1:length(chirp_phases);
                [p, s, mu] = polyfit(t, chirp_phases, 9);
                trend = polyval(p, t, [], mu);
                corrected_phases = chirp_phases - trend;
                frame(:,kk,rangeBinIndex) = abs(chirps_at_rx_rangebin).*exp(1j*corrected_phases);
                
%                 figure;
%                 subplot(2,1,1); plot(chirp_phases); hold on; plot(trend);
%                 subplot(2,1,2); plot(corrected_phases);
            end        
            radarCube.data{ii} = frame;
        end
    end
    
    newradarCube = radarCube;
end