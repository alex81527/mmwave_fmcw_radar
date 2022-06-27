function [out] = decompose_sincs(adc_data, params, numSincs, ...
    sincWindowSize, out_given, plotting)
    if ~isempty(out_given)
        out = out_given;
        return
    end

    fftsize = 512; %params.opRangeFFTSize;
    w = 2*pi*[0:fftsize-1]/fftsize;
    rangefft_output = fft(adc_data, fftsize);
    if iscolumn(rangefft_output)
        rangefft_output = rangefft_output.';
    end
    n1 = norm(rangefft_output);
    
    cf=0;
    out.sincs = zeros(numSincs, fftsize);
%     out.search_sincs = zeros(fftsize, 2*params.w_steps+1, numSincs);
    out.rangefft_idx = zeros(1, numSincs);
    out.amps = zeros(1, numSincs);
    out.phis = zeros(1, numSincs);
    out.dists = zeros(1, numSincs);
    for ii=1:numSincs
        [M, I] = max(abs(rangefft_output));
        ct = w(I);
        phi = angle(rangefft_output(I));
        peak_sinc = M*exp(1j*phi)*exp(-1j*(w-ct)*(sincWindowSize-1)/2).*diric(w - ct, sincWindowSize);
        out.sincs(ii, :) = peak_sinc;
        out.rangefft_idx(ii) = I;
        out.amps(ii) = M;
        out.phis(ii) = phi;
        out.dists(ii) = (ct*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
        
%         ct_idx = (I-1)*params.w_steps +1 ;
%         ct = params.w(ct_idx);
%         cand_ct_idx = [-params.w_steps:params.w_steps] + ct_idx;
%         cand_ct = params.w(cand_ct_idx);
%         out.search_sincs(:,:,ii) = getChirpSincs(params, cand_ct);
 
        cf = cf + 1;
        if plotting
            figure(cf);
            subplot(2,1,1); plot(abs(rangefft_output)); hold on; plot(abs(peak_sinc)); %xlim([0 1000]);
            subplot(2,1,2); plot(abs(rangefft_output-peak_sinc)); %xlim([0 1000]);    
            n2 = norm(rangefft_output-peak_sinc);
            fprintf('[%d] I=%d %.1f/%.1f (%.1f%%)\n', cf, I, n2, n1, 100*n2/n1);
        end
        rangefft_output = rangefft_output - peak_sinc;
    end
%     out.dists = (w(out.rangefft_idx)*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
end