function f = obj_func(x)
    mid = (length(x)-1)/3;
    amps = x(1:mid);
    dists = x(mid+1:2*mid);
    phis = x(2*mid+1:3*mid);
    const_phase = x(end);
    
    freqSlope = 64.9850; %params.freqSlope;
    sampleRate = 18.5; %params.sampleRate;
    numSamplePerChirp = 1024; %params.numSamplePerChirp;
    fftsize = 1024*8; %length(t);
%     amps = [0.5486 3.891 0.8065 0.632];%[0.8 2.8 0.84 0.92 0.44];
%     dists = [0.079 0.7032 0.8928 1.362];%[0.072 0.5942 1.11 1.215 1.314];%[150.5 ]*params.rangeResolutionsInMeters;
    IF_freqs = freqSlope*1e12*2/3e8*dists;
%     delta_f = 256*freqSlope/sampleRate*1e6;
%     phis = 4*pi*dists/(3e8/(77e9+0*delta_f));
%     t = [0:numSamplePerChirp-1].'/(sampleRate*1e6);
%     sinusoids = amps.*exp(1j*2*pi.*IF_freqs.*t).*exp(1j*phis).*exp(1j*const_phase);
%     sig_t = sum(sinusoids, 2) + 0.0*(randn(length(t),1) + 1j*randn(length(t),1));
    % figure;
%     sig_f = fft(sig_t, fftsize).';
    load('adc_data.mat');
    rangefft_output = fft(double(adc_data), fftsize).'; 
    % stem(abs(sig_f));
    w = 2*pi*[0:fftsize-1]/fftsize;
    center = IF_freqs/(sampleRate*1e6)*2*pi;
    D = zeros(1, length(w));
    for ii = 1:length(center)
        ct = center(ii);
        window_size = length(adc_data);
        D = D + amps(ii)*exp(1j*phis(ii))*exp(-1j*(w-ct)*(window_size-1)/2).*diric(w - ct, window_size);
        [M, I] = min(abs(w-ct));
        peak_idx(ii) = I;
    end
%     figure;
%     subplot(2,1,1); plot(abs(rangefft_output)); hold on; plot(abs(D)); xlim([0 600]);
%     subplot(2,1,2); plot(angle(rangefft_output)); hold on; plot(angle(D)); xlim([0 600]);
%     figure;
%     plot(angle(exp(1j*phis))); hold on; plot(angle(sig_f(peak_idx)), 'ro');
    f = norm(rangefft_output - D)^2;
end