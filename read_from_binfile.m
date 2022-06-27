function [adcRawOutput] = read_from_binfile(numBinFile, nFrame, numChirpPerFrame, numSamplePerChirp)
    bytePerSample = 4;
    lanes = 4;
%     numSamplePerChirp = 1024;
%     numChirpPerFrame = 16*128;
    frameSize = numChirpPerFrame*numSamplePerChirp*lanes*bytePerSample;
    adcRawOutput = zeros(numSamplePerChirp, numChirpPerFrame, nFrame);
    byteNeeded = 0;
    frameIdx =0;
    for ii=0:numBinFile-1
        filename = sprintf('C:/ti/mmwave_studio_02_01_01_00/mmWaveStudio/PostProc/adc_data_Raw_%d.bin', ii);
        [fid, errmsg] = fopen(filename, 'r');
        ff = dir(filename);
        totalbytes = ff.bytes;
        
        byteIdx = 0;
        if byteNeeded~=0
            rawdata = fread(fid, byteNeeded/2,  'uint16=>single');
            rawdata = [last_rawdata ; rawdata];
            byteIdx = byteIdx + byteNeeded;
            byteNeeded = 0;
            rawdata = rawdata - ( rawdata >=2.^15).* 2.^16;
            rawData8 = reshape(rawdata, [8, length(rawdata)/8]);
            rawDataI = reshape(rawData8(2,:), [], 1);
            rawDataQ = reshape(rawData8(6,:), [], 1);
            frameData = rawDataI + 1j*rawDataQ;
            frameIdx = frameIdx +1;
            adcRawOutput(:, :, frameIdx) = reshape(frameData, numSamplePerChirp, numChirpPerFrame); 
        end
        
        while byteIdx < totalbytes
            if  byteIdx + frameSize <= totalbytes
                rawdata = fread(fid, frameSize/2,  'uint16=>single');  % reads 2 bytes each time
                byteIdx = byteIdx + frameSize;
            else
                readbyte = totalbytes - byteIdx;
                byteNeeded = frameSize - readbyte;
                last_rawdata = fread(fid, readbyte/2,  'uint16=>single');
                if ii~=numBinFile-1
                    break;
                else
                    error('Frame incomplete.');
                end
            end
            % signed value adjustment
            rawdata = rawdata - ( rawdata >=2.^15).* 2.^16;
            % Convert 4 lane LVDS data to one matrix
            rawData8 = reshape(rawdata, [8, length(rawdata)/8]);
%             rawDataI = reshape(rawData8(1:4,:), [], 1);
%             rawDataQ = reshape(rawData8(5:8,:), [], 1);
            rawDataI = reshape(rawData8(2,:), [], 1);
            rawDataQ = reshape(rawData8(6,:), [], 1);
            frameData = rawDataI + 1j*rawDataQ;
            frameIdx = frameIdx +1;
            adcRawOutput(:, :, frameIdx) = reshape(frameData, numSamplePerChirp, numChirpPerFrame); 
        end
        fclose(fid);
    end
    assert(frameIdx == nFrame);
end