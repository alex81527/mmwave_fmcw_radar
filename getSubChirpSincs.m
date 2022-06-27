function [out] = getSubChirpSincs(params, ct)
    if params.subChirpSincHashMap.isKey(ct)
        out = params.subChirpSincHashMap(ct);
    else
        idx = [1:params.w_steps:length(params.w)]-1;
        shift = floor(ct/params.w_stepsize);
        idx = mod(idx-shift,length(params.w))+1;
        params.subChirpSincHashMap(ct) = params.subChirpSincs(idx);
        out = params.subChirpSincHashMap(ct);
    end
end