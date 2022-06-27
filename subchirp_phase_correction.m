function [newradarCube] = subchirp_phase_correction(radarCube, rangeBinIndex, correction_mode)
    N_frames = 1; %length(radarCube.data);
    N_chirps = radarCube.dim.numChirps;
    N_rxchains = 1; %radarCube.dim.numRxChan;
    
    stride = radarCube.dim.subChirpStrideStep;
    delta_f = stride*(1/radarCube.rfParams.sampleRate)*radarCube.rfParams.freqSlope*1e6;
    object_distance = (rangeBinIndex-1)*radarCube.rfParams.rangeResolutionsInMeters;
    theoretical_phase_slope = 4*pi*(object_distance)/3e8*delta_f;

    if strcmp(correction_mode, 'theoretical')
        for ii = 1:N_frames
            frame = radarCube.data{ii};
            for jj = 1:N_chirps
                for kk = 1:N_rxchains
                    subchirps_at_rx_rangebin = squeeze(frame(jj,:,kk,rangeBinIndex));
                    theoretical_correction = exp(-1j*theoretical_phase_slope*[0:radarCube.dim.numSubChirps-1]);
                    corrected_subchirps = subchirps_at_rx_rangebin .* theoretical_correction;                                   
                end
                frame(jj,:,kk,rangeBinIndex) = corrected_subchirps;
            end
            radarCube.data{ii} = frame;
        end
    elseif strcmp(correction_mode, 'fitted')
        for ii = 1:N_frames
            frame = radarCube.data{ii};
            for jj = 1:N_chirps
                for kk = 1:N_rxchains
                    subchirps_at_rx_rangebin = squeeze(frame(jj,:,kk,rangeBinIndex));
                    unwrap_phases = unwrap(angle(subchirps_at_rx_rangebin));
                    t = 1:length(unwrap_phases);
                    [p, s, mu] = polyfit(t, unwrap_phases, 1); % linear fit
%                     fitted_points = polyval(p, t, [], mu);
                    fitted_slope = p(1)/mu(2);
                    fitted_correction = exp(-1j*fitted_slope*[0:radarCube.dim.numSubChirps-1]);
                    corrected_subchirps = subchirps_at_rx_rangebin .* fitted_correction;                                   
                end
                frame(jj,:,kk,rangeBinIndex) = corrected_subchirps;
            end
            radarCube.data{ii} = frame;
        end
    elseif strcmp(correction_mode, 'characterization')
        frame = radarCube.data{1}; % frame 1
        subchirps_at_rx_rangebin = squeeze(frame(1,:,1,rangeBinIndex)); % chirp 1, rx 1
        unwrap_phases = unwrap(angle(subchirps_at_rx_rangebin));
        t = 1:length(unwrap_phases);
        [p, s, mu] = polyfit(t, unwrap_phases, 1); % linear fit
        fitted_points = polyval(p, t, [], mu);
        fitted_slope = p(1)/mu(2);
        fitted_correction = exp(-1j*fitted_slope*[0:radarCube.dim.numSubChirps-1]);
        theoretical_correction = exp(-1j*theoretical_phase_slope*[0:radarCube.dim.numSubChirps-1]);
        figure;
        plot(unwrap_phases, 'b'); hold on;
        plot(p(1)*(t-mu(1))/mu(2)+p(2), 'ro'); hold on;
        plot(unwrap_phases(1) + theoretical_phase_slope*[0:radarCube.dim.numSubChirps-1], 'k');
        legend('measured', sprintf('fitted slope %.5f', fitted_slope), ...
            sprintf('theoretical slope %.5f', theoretical_phase_slope), ...
            'Location', 'southeast', 'FontSize', 12);
        title(sprintf('phases across subchirps'));
        fprintf('theoretical slope %.5f, fitted slope %.5f, dev: %.3f %%\n',theoretical_phase_slope,...
            fitted_slope, 100*(fitted_slope-theoretical_phase_slope)/theoretical_phase_slope);
        
        
        first_subchirp_phase = mod(angle(subchirps_at_rx_rangebin(1)), 2*pi);
        theo_corr_phases = mod(angle(subchirps_at_rx_rangebin.* theoretical_correction), 2*pi);
        fit_corr_phases = mod(angle(subchirps_at_rx_rangebin.* fitted_correction), 2*pi);
        
        figure;
        subplot(3,1,1); plot(angle(subchirps_at_rx_rangebin)); title('original phases');
        subplot(3,1,2); plot(theo_corr_phases); hold on; 
        plot([1 length(subchirps_at_rx_rangebin)], [first_subchirp_phase first_subchirp_phase]);
        title('theoretical correction'); legend('ideal', 'correction result');
        subplot(3,1,3); plot(fit_corr_phases); hold on; 
        plot([1 length(subchirps_at_rx_rangebin)], [first_subchirp_phase first_subchirp_phase]);
        title('fitted correction'); legend('ideal', 'correction result');
        
        fprintf('Mean deviation: theoretical %.3f, fitted %.3f\n', ...
            mean(abs(theo_corr_phases-first_subchirp_phase)), ...
            mean(abs(fit_corr_phases-first_subchirp_phase)));
    else
        error('NotImplemented');
    end
    
    newradarCube = radarCube;
end