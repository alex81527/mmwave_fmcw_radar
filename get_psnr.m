function [out] = get_psnr(audio_samples, fs)
    wl = fix(0.4*fs);
    nov = fix(wl/2);
%     y2 = reshape(audio_samples(1:end-mod(length(audio_samples),wl)), wl, []);
    y2 = buffer(audio_samples, wl, nov);
    y3 = mean(abs(y2).^2, 1);
    out = 10*log10(max(y3)/min(y3));
%     nff = 2^nextpow2(wl);
%     yyy = fft(reshape(audio_samples(1:end-mod(length(audio_samples),wl)), wl, []),nff);
%     zzz = max(abs(yyy), [],1);
    
end