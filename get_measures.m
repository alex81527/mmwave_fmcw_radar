function [llr, wss, stoi_val, nist_snr, snr_seg] = get_measures(y_clean, fs1, y_enhanced, fs2, do_plot)
    addpath('./snreval-master');
    target_samprate = 16000;
    if fs1 > target_samprate
        y_clean = resample(y_clean, target_samprate, fs1);
    end
    if fs2 > target_samprate
        y_enhanced = resample(y_enhanced, target_samprate, fs2);
    end
    x           = y_clean(:);                             % clean speech column vector
    y           = y_enhanced(:);                      % processed speech column vector
    
    % clean speech is shorter since y_enhanced includes leading noise samples 
    % find out how many zeros needed to prepend the clean speech to align
    % them
    [x2, y2] = alignsignals(x,y);
    zeros_prepended = length(x2) - length(x);
    x2= x2(  zeros_prepended + [1: min(length(y)-zeros_prepended,length(x))] );
    y2= y2(  zeros_prepended + [1: min(length(y)-zeros_prepended,length(x))] );
%      length(x)/target_samprate
    if do_plot
        figure; 
        subplot(2,1,1); plot(xcorr(x,y));
    %     subplot(3,1,1); plot([1:length(x3)]/Fs2, x3);
    %     subplot(3,1,2); plot([1:length(x2)]/Fs2, x2);
        subplot(2,1,2); plot(x2./max(abs(x2))); hold on; plot(y2./max(abs(y2)));
    end
    % get_psnr(x2,Fs2)
    % get_psnr(x3, Fs2)
    [snr_all, snr_seg] = comp_snr(x2,y2, target_samprate);
    llr = comp_llr(x2,y2, target_samprate);
    wss =comp_wss(x2,y2, target_samprate);
    stoi_val = stoi(x2,y2, target_samprate);
    nist_snr = nist_stnr_m(y2,target_samprate);
    fprintf('LLR=%.2f, WSS=%.2f, Intelligibility score=%.2f, NIST-SNR=%.2f\n',llr, wss, stoi_val,nist_snr);
    fprintf('snr_all=%.2f, snr_seg=%.2f\n',snr_all, snr_seg);
end