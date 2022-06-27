function [out] = get_ssnr(audio_samples, fs)
    wl = fix(0.3*fs);
    y2 = reshape(audio_samples(1:end-mod(length(audio_samples),wl)), wl, []);
    y3 = mean(abs(y2).^2, 1);
    out = 10*log10(mean(y3(y3>max(y3)*0.1)) / mean(mink(y3,2)));
%     figure; plot(y3);
% figure; plot(audio_samples);
%     nff = 2^nextpow2(wl);
%     yyy = fft(reshape(audio_samples(1:end-mod(length(audio_samples),wl)), wl, []),nff);
%     zzz = max(abs(yyy), [],1);
    
end