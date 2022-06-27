function [params] = read_from_json(filename)
    [fid, errmsg] = fopen(filename, 'r');    
    a  = fread(fid);
    b = jsondecode(char(a.'));
    chirpProfileId = b.mmWaveDevices.rfConfig.rlChirps.rlChirpCfg_t.profileId;
    chirpProfile = b.mmWaveDevices.rfConfig.rlProfiles(chirpProfileId+1).rlProfileCfg_t;
   
    params.startFreq = chirpProfile.startFreqConst_GHz*1e9;
    params.freqSlope = chirpProfile.freqSlopeConst_MHz_usec*1e12;
    params.chirpCycleTime = (chirpProfile.idleTimeConst_usec + chirpProfile.rampEndTime_usec)*1e-6;
    params.sampleRate = chirpProfile.digOutSampleRate*1e3;
    params.numSamplePerChirp = chirpProfile.numAdcSamples;
    params.bandwidth = params.numSamplePerChirp*params.freqSlope/params.sampleRate;
    params.rangeResolutionsInMeters = 3e8/2/params.bandwidth;
%     params.dopplerResolutionMps = radarCube.rfParams.dopplerResolutionMps;
    params.framePeriodicity = b.mmWaveDevices.rfConfig.rlFrameCfg_t.framePeriodicity_msec*1e-3;
%     params.numFrames = radarCube.dim.numFrames;
    params.numChirps = (b.mmWaveDevices.rfConfig.rlFrameCfg_t.chirpEndIdx+1)*...
        b.mmWaveDevices.rfConfig.rlFrameCfg_t.numLoops;
%     params.numRxChan = radarCube.dim.numRxChan;
    
    params.opRangeFFTSize = 2^13; %max(2^13, 2^nextpow2(params.numSamplePerChirp));
    params.opDopplerFFTSize = 2^13; %max(2^13, 2^nextpow2(params.numChirps));
    
    fclose(fid);
end