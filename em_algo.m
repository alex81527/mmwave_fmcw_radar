function [out2] = em_algo(adc_data, params, numSincs, sincWindowSize,...
    out1, plotting)
    empty_flag = isempty(out1);
    % initialize 
    fftsize = params.opRangeFFTSize;
    w = 2*pi*[0:fftsize-1]/fftsize;
%     if empty_flag
%         [out1] = decompose_sincs(adc_data, params, numSincs,...
%         sincWindowSize, [], plotting);
% %         out1.sincs = zeros(numSincs, fftsize);
%     end
    
    rangefft_output = fft(adc_data, fftsize);
    rangefftMaxAbs = max(abs(rangefft_output));
    threshold = 0.1*rangefftMaxAbs;
    if iscolumn(rangefft_output)
        rangefft_output = rangefft_output.';
    end
    
    if empty_flag
        out1.sincs = zeros(numSincs, fftsize);
    end
    leftover = rangefft_output - sum(out1.sincs, 1);
    sincs_copy = out1.sincs; 
    init_norm = norm(leftover);
    maxIter = 50;
    out2.sincs = zeros(numSincs, fftsize);
    out2.rangefft_idx = zeros(1, numSincs);
    out2.w_idx = zeros(1, numSincs);
    out2.amps = zeros(1, numSincs);
    out2.phis = zeros(1, numSincs);
    out2.dists = zeros(1, numSincs);
    round_amps = zeros(maxIter, numSincs);
    round_w_idx = zeros(maxIter, numSincs);
    round_phis = zeros(maxIter, numSincs);
    round_norm = zeros(1, maxIter);
    for kk=1:maxIter
        for ii=1:numSincs
            sig = sincs_copy(ii,:) + leftover;
            
%             fprintf('%.1f %.1f\n', max(abs(sig)), threshold);
%             if max(abs(sig)) < threshold
%                 break;
%             end
            
            if empty_flag 
                if kk==1           
                    [M, I] = max(abs(sig));
                    ct_idx = (I-1)*params.w_steps +1 ;
                    ct = params.w(ct_idx);
                    out2.sincs(ii, :) = sig(I).*getChirpSincs(params, ct);
                    out2.rangefft_idx(ii) = I;
                    out2.w_idx(ii) = ct;
                    out2.amps(ii) = abs(sig(I));
                    out2.phis(ii) = angle(sig(I));
                    out2.dists(ii) = (ct*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
                else
                    [M, I] = max(abs(sig));
                    ct_idx = (I-1)*params.w_steps +1 ;
                    ct = params.w(ct_idx);
                    cand_ct_idx = [-params.w_steps:params.w_steps] + ct_idx;
                    cand_ct = params.w(cand_ct_idx);
                    w_diff = ct - cand_ct;

                    adjust = 1./...
                        (exp(-1j*(w_diff)*(sincWindowSize-1)/2).*...
                        sin(sincWindowSize*w_diff/2)./(sincWindowSize*sin(w_diff/2)));

                    search_sincs = getSearchChirpSincs(params, ct);
                    errors = sum(abs(sig.' - sig(I).*adjust.*search_sincs),1);

    %                 figure; plot(errors);
                    [M2,I2] = min(errors);

                    peak_gain = sig(I)*adjust(I2);
                    ct = cand_ct(I2);

                    out2.sincs(ii, :) = peak_gain.*search_sincs(:,I2).'; 
                    out2.rangefft_idx(ii) = I;
                    out2.w_idx(ii) = cand_ct_idx(I2);
                    out2.amps(ii) = abs(peak_gain);
                    out2.phis(ii) = angle(peak_gain);
                    out2.dists(ii) = (ct*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
                end
            else
%                 [M, I] = max(abs(sig));
                I = out1.rangefft_idx(ii);
                ct_idx = (I-1)*params.w_steps +1 ;
                ct = params.w(ct_idx);
                known_ct = params.w(out1.w_idx(ii));
                w_diff = ct - known_ct;
                adjust = 1./...
                    (exp(-1j*(w_diff)*(sincWindowSize-1)/2).*...
                    sin(sincWindowSize*w_diff/2)./(sincWindowSize*sin(w_diff/2)));
                peak_gain = sig(I)*adjust;
                
                if sincWindowSize==params.chirpSincWindowSize
%                     error('!!');
                    out2.sincs(ii, :) = peak_gain.*getChirpSincs(params, known_ct); 
                else
                    out2.sincs(ii, :) = peak_gain.*getSubChirpSincs(params, known_ct); 
                end
                    
                out2.rangefft_idx(ii) = I;
                out2.w_idx(ii) = out1.w_idx(ii);
                out2.amps(ii) = abs(peak_gain);
                out2.phis(ii) = angle(peak_gain);
                out2.dists(ii) = (known_ct*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
                       
            end

            leftover = leftover + sincs_copy(ii,:) - out2.sincs(ii, :);
%             
%             figure;
%             plot(abs(sig)); hold on;
%             plot(abs(out2.sincs(ii,:))); hold on;
%             plot(abs(sig-out2.sincs(ii,:)));
%             legend('sig', 'sincs', 'sig - sincs');
%             plot([1 fftsize], threshold*ones(1,2));
% 
%             xlim([max(I-200,0) I+200])
%             fprintf('[iter=%d][sinc%d][rangebin=%d] leftover norm: %.1f\n', ...
%                 kk,ii,out2.rangefft_idx(ii),norm(leftover));
        end
        sincs_copy = out2.sincs;
        round_amps(kk, :) = out2.amps;
        round_w_idx(kk, :) = out2.w_idx;
        round_phis(kk, :) = out2.phis;
        round_norm(kk) = norm(rangefft_output - sum(out2.sincs, 1));
        
        if round_norm(kk) > init_norm
%             error('round_norm(kk) > init_norm');
            fprintf('round_norm(kk) > init_norm, ends with %d rounds\n', kk);
            break;
        end
        
        if kk>1 && 100*(round_norm(kk-1) - round_norm(kk))/init_norm < 0.1
            break
        end
    end
    
    if plotting
        figure;
        subplot(4,1,1);
        plot(100*[init_norm round_norm(1:kk)]./init_norm); 
        title('leftover norm percentage');
        subplot(4,1,2);
        plot(round_amps(2:kk,:) - round_amps(1,:)); title('amps diff');
        subplot(4,1,3);
        plot(round_w_idx(2:kk,:) - round_w_idx(1,:)); title('w idx diff');
        subplot(4,1,4);
        plot(round_phis(2:kk,:) - round_phis(1,:)); title('phis diff');
        
        figure;
        subplot(2,1,1); plot(abs(rangefft_output), 'k'); hold on; 
        plot(abs(sum(out1.sincs,1)),'b'); hold on; 
        plot(abs(sum(out2.sincs,1)),'r'); %xlim([0 1000]);
        legend('rangefft', 'initial', 'refined');
        subplot(2,1,2); plot(angle(rangefft_output), 'k'); hold on; 
        plot(angle(sum(out1.sincs,1)),'b'); hold on; 
        plot(angle(sum(out2.sincs,1)),'r'); %xlim([0 1000]);
        legend('rangefft', 'initial', 'refined');
        
        figure;
        plot(abs(rangefft_output), 'k'); hold on;
        plot(abs(rangefft_output - sum(out2.sincs,1)), 'r'); %xlim([0 1000]);
        legend('rangefft', 'leftover');
%         plot(1:fftsize, abs(new_sincs)); xlim([0 1000]);
    end
end