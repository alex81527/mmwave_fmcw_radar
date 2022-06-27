function [fft_output, I] = rangeFFT(adc_data, params)
    fft_output = fft(adc_data, params.opRangeFFTSize);
    freq_interval = params.sampleRate/params.opRangeFFTSize;
    dist_interval = freq_interval*3e8/(2*params.freqSlope);
    % find objects within 0.1 - 1 m
    min_distance = 0.1;
    max_distance = 1.0;
    min_idx = ceil(min_distance/dist_interval);
    max_idx =  ceil(max_distance/dist_interval);
    [M, I] = max(abs(fft_output(min_idx:max_idx)));
    % [M, I] = max(abs(rangefft_output(1:end)));
    I = I + (min_idx-1) ;
    object_dist = dist_interval * (I-1);
    fprintf('Within %.2f - %.2f m, object found at %.3f m index=%d \n', min_distance, max_distance,...
        object_dist, I);
    figure;
    subplot(3,1,1); plot(dist_interval*[0:params.opRangeFFTSize-1], abs(fft_output));
    xlabel('range (m)'); ylabel('range FFT output (dB)'); title('Range FFT'); xlim([min_distance max_distance]);
    subplot(3,1,2); plot(db(fft_output));
    xlabel('index'); ylabel('range FFT output (dB)');
    subplot(3,1,3); plot(abs(fft_output)); hold on;
    xlabel('index'); ylabel('range FFT output (abs)'); xlim([1 max_idx]);
    stem(I, M);
end