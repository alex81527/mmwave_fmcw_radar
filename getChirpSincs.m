function [out] = getChirpSincs(params, ct)
    if params.chirpSincHashMap.isKey(ct)
        out = params.chirpSincHashMap(ct);
    else
        idx = [1:params.w_steps:length(params.w)]-1;
        shift = floor(ct/params.w_stepsize);
        idx = mod(idx-shift,length(params.w))+1;
        params.chirpSincHashMap(ct) = params.chirpSincs(idx);
        out = params.chirpSincHashMap(ct);
    end
end