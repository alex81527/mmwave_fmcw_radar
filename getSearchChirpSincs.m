function [out] = getSearchChirpSincs(params, ct)
    if params.chirpSearchSincHashMap.isKey(ct)
        out = params.chirpSearchSincHashMap(ct);
    else
        cand_ct = ct + params.w_stepsize*[-params.w_steps:params.w_steps];
        idx = [1:params.w_steps:length(params.w)].'-1;
        shift = floor(cand_ct/params.w_stepsize);
        idx = mod(idx-shift,length(params.w))+1;
        params.chirpSearchSincHashMap(ct) = params.chirpSincs(idx);
        out = params.chirpSearchSincHashMap(ct);
    end
end