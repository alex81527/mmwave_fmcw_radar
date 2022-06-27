%%
amps = [0.8 2.8 0.84 0.92 0.44];%[0.5486 3.891 0.8065 0.632]*1e6;
dists = [0.072 0.5942 1.11 1.215 1.314];%[0.079 0.7032 0.8928 1.362];%[150.5 ]*params.rangeResolutionsInMeters;
IF_freqs = params.freqSlope*1e12*2/3e8*dists;
delta_f = 256*params.freqSlope/params.sampleRate*1e6;
phis = 4*pi*dists/(3e8/(77e9+0*delta_f));
t = [0:params.numSamplePerChirp-1].'/(params.sampleRate*1e6);
fftsize = 1024*8; %length(t);
sinusoids = amps.*exp(1j*2*pi.*IF_freqs.*t).*exp(1j*phis);
sig_t = sum(sinusoids, 2) + 0.0*(randn(length(t),1) + 1j*randn(length(t),1));
% figure;
sig_f = fft(sig_t, fftsize).';
% stem(abs(sig_f));
w = 2*pi*[0:fftsize-1]/fftsize;
center = IF_freqs/(params.sampleRate*1e6)*2*pi;
D = zeros(1, length(w));
for ii = 1:length(center)
    ct = center(ii);
    window_size = length(sig_t);
    D = D + amps(ii)*exp(1j*phis(ii))*exp(-1j*(w-ct)*(window_size-1)/2).*diric(w - ct, window_size);
    [M, I] = min(abs(w-ct));
    peak_idx(ii) = I;
end

% figre;
% plot(abs(rangefft_output)); hold on; plot(abs(D)); xlim([0 500]);
% plot(abs(fft(adc_data(1:window_size).', fftsize))); hold on; plot(window_size*abs(D)); xlim([0 500]);
% plot(abs(sig_f)); hold on; plot(length(sig_t)*abs(D)); xlim([0 500]);

%%
close all;
load('adc_data.mat');
fftsize = params.opRangeFFTSize;
numSincs = 15;
scaling = 1;
diricWindowSize = length(adc_data);
w = 2*pi*[0:fftsize-1]/fftsize;
rangefft_output = fft(double(adc_data), fftsize).';
[amps, phis, pk_idx, sincs] = decompose_sincs(adc_data, fftsize, numSincs, diricWindowSize, scaling, 1);
d = (w(pk_idx)*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
figure; plot(abs(rangefft_output)); hold on; plot(abs(sum(sincs,1))); xlim([0 2500]);

angle(sincs(:, pk_idx(3)))
angle(sum(sincs(:,pk_idx(3)),1))

abs(sincs(:, pk_idx(3)))
abs(sum(sincs(:,pk_idx(3)),1))

figure;
plot(angle(sum(sincs(:,pk_idx),1))); hold on; 
plot(angle(rangefft_output(pk_idx)), 'ro'); hold on;
plot(angle(diag(sincs(:,pk_idx))), 'k'); hold on; 
figure;
sorted_idx = sort(pk_idx);
sorted_d = (w(sorted_idx)*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
plot(unwrap(angle(exp(1j*4*pi*sorted_d/(3e8/77e9)))), 'm'); hold on;
plot(unwrap(angle(rangefft_output(sorted_idx))), 'bo');


% diricWindowSize256 = 256;
% rangefft_output256 = fft(adc_data(1:diricWindowSize256), fftsize);
% centers = w(pk_idx).';
% sinc256 = (1/4)*(amps.*exp(1j*phis)).'.*exp(-1j*(w-centers)*(diricWindowSize256-1)/2).*diric(w - centers, diricWindowSize256);
% figure; plot(abs(rangefft_output256)); hold on; plot(abs(sum(sinc256,1)));
% 
% angle(diag(sinc256(:, pk_idx))).'
% angle(sum(sinc256(:,pk_idx),1))
% figure;
% plot(angle(sum(sinc256(:,pk_idx),1))); hold on; 
% plot(angle(rangefft_output256(pk_idx)), 'ro'); hold on;
% plot(angle(diag(sinc256(:,pk_idx))), 'k'); hold on; 
% d = (w(pk_idx)*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
% % plot(angle(exp(1j*4*pi*d/(3e8/77e9))), 'm'); 

clear sinc256 rangefft_output256
for ii=1:4
    diricWindowSize256 = 128;
    centers = w(pk_idx).';
    phase_adjust = 4*pi*d/3e8*((ii-1)*diricWindowSize256*params.freqSlope/params.sampleRate*1e6);
    sinc256(:,:,ii) = (1/4)*(amps.*exp(1j*phis).*exp(1j*phase_adjust)).'.*exp(-1j*(w-centers)*(diricWindowSize256-1)/2).*diric(w - centers, diricWindowSize256);
    data_in = adc_data((ii-1)*diricWindowSize256+[1:diricWindowSize256]);
%     [amps2, phis2, pk_idx2, sincs2] = decompose_sincs(data_in, fftsize, numSincs, 384, scaling, 0);
    rangefft_output256(ii, :) = fft(data_in, fftsize);
%     d2 = (w(pk_idx2)*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
    
%     subchirp2_amps(ii,:) = amps2; 
%     subchirp2_phis(ii,:) = phis2; 
%     subchirp2_pk_idx(ii,:) = pk_idx2; 
%     subchirp2_d(ii,:) = d2; 
    
    sum_sincs = sum(sinc256(:,:,ii),1);
    
    figure(55); 
    subplot(4,1,ii);
    plot(abs(rangefft_output256(ii, :)), 'b'); hold on; 
    stem(pk_idx, amps/4, 'k'); hold on;
    plot(abs(sum_sincs), 'r'); xlim([0 2500]);
    
    figure(56); 
    subplot(4,1,ii);
    plot(angle(rangefft_output256(ii, :)), 'b'); hold on; 
    plot(pk_idx, angle(sum_sincs(pk_idx)), 'ko'); hold on;
    plot(angle(sum_sincs), 'r'); xlim([0 2500]);
end

verify_idx = 3;
for ii=1:4
    sum_sincs = sum(sinc256(:,:,ii),1);
    
    figure(65);
    subplot(4,1,ii);
    plot(abs(sinc256(:, pk_idx(verify_idx),ii))); hold on;
    plot(verify_idx, abs(sinc256(verify_idx, pk_idx(verify_idx),ii)), 'k^');
    plot(verify_idx, abs(rangefft_output256(ii, pk_idx(verify_idx))), 'ro');
    title(sprintf('subchirp sincs abs at %.3f m', d(verify_idx)));

    figure(66);
    subplot(4,1,ii);
    plot(angle(sinc256(:, pk_idx(verify_idx),ii))); hold on;
    plot(verify_idx, angle(sinc256(verify_idx, pk_idx(verify_idx),ii)), 'k^');
    plot(verify_idx, angle(rangefft_output256(ii, pk_idx(verify_idx))), 'ro');
    title(sprintf('subchirp sincs angle at %.3f m', d(verify_idx)));
end

subchirp2_pk_idx
subchirp2_d
subchirp2_phis(end-6:end,:)
4*pi*subchirp2_d(end,:)/3e8*(params.freqSlope/params.sampleRate*1e6)


for ii=1:4
    [amps3, phis3, pk_idx3, sincs3] = decompose_sincs(adc_data((ii-1)*256+[1:256]), fftsize, numSincs, 256, scaling, 0);
    rangefft_output256 = fft(adc_data((ii-1)*256+[1:256]), fftsize);
    d3 = (w(pk_idx3)*(params.sampleRate*1e6)/2/pi)*3e8/2/(params.freqSlope*1e12);
    
    subchirp3_amps(ii,:) = amps3; 
    subchirp3_phis(ii,:) = phis3; 
    subchirp3_pk_idx(ii,:) = pk_idx3; 
    subchirp3_d(ii,:) = d3; 
    figure(66); 
    subplot(4,1,ii);
    plot(abs(rangefft_output256)); hold on; plot(abs(sum(sincs3,1))); xlim([0 2500]);
end
sort(d)
sort(subchirp3_d, 2)
% b = rangefft_output(cand_peak_idx).';
% A = cand_sincs(:, cand_peak_idx).';
% x = inv(A)*b;
% abs(x)
% angle(x)
% norm(sig_f - sum(cand_sincs,1))
% norm(sig_f - sum(cand_sincs.*x,1))
% figure;
% plot(abs(sig_f)); hold on; plot(abs(sum(cand_sincs,1))); hold on; plot(abs(sum(cand_sincs.*x,1)));


% [pks,locs,w,p] = findpeaks(abs(sig_f), 'MinPeakHeight', 400);
% figure;
% subplot(2,1,1); plot(abs(sig_f)); hold on; plot(length(sig_t)*abs(D)); xlim([0 200]);
% subplot(2,1,2); plot(angle(sig_f)); hold on; plot(angle(D));
% [M, I] = max(abs(sig_f - length(sig_t)*D))

%%
N = 1024;
n = 0:N-1;

w0 = 2*pi/5;
x = sin(w0*n)+10*sin(2*w0*n);
s = spectrogram(x);

spectrogram(x,'yaxis')

%%
% [y, Fs] = audioread('Original.wav');
% [y, Fs] = audioread('IKEABag.wav');

figure;
for ii=1:9
    [y, Fs] = audioread(sprintf('myvoice%d.wav', ii));
    % sound(y, Fs);
    Nx = length(y);
    nsc = floor(Nx/16);
    nov = floor(nsc/2);
    nff = max(256,2^nextpow2(nsc));
    subplot(3,3,ii);
    spectrogram(y,hamming(nsc),nov,nff, Fs, 'yaxis'); %,'MinThreshold',-70)
    ylim([0 0.9]); title(sprintf('Digit %d', ii));
end
% ylim([0 0.25])
% spectrogram(y,'yaxis')
% figure; plot(y);
% figure; plot(abs(fft(y(23000+[1:1024]))));
%%
fs = 1e3;
t = 0:1/fs:1;
x = [1 2]*sin(2*pi*[50 250]'.*t) + randn(size(t))/10;
highpass(x,150,fs)

fs = 1/44e-6;
t = [0:255]./fs;
resolution = 1/(255/fs)
x = cos(2*pi*440.*t) + 0.*rand(1, 256);
x2 = cos(2*pi*440.*t+0.5*pi) + 0.*rand(1, 256);
X = fft(x); angle(X(6))
X2 = fft(x2); angle(X2(6))
delta = angle(X(6)) - angle(X2(6));
X2(6) = X2(6)*exp(1j*delta);
X2(252) = X2(252)*exp(1j*delta);
x2_prime = ifft(X2);

figure; plot(x); hold on; plot(x2); hold on; plot(real(x2_prime));
figure; plot(abs(X));

%%
fs = 5e3;
t = 0:1/fs:1;
x = [2 1 2]*sin(2*pi*[50 150 250]'.*t) + randn(size(t))/10;
filtered_x = bandpass(x,[100 200],fs);
figure;
plot(x); hold on; plot(filtered_x);

figure;
plot(fs/1024*[0:1023],db(fft(x, 1024))); hold on; 
plot(fs/1024*[0:1023],db(fft(filtered_x, 1024)));

%%
b = 1400/5*67e-6*2048;
bpFilt = designfilt('bandpassiir','FilterOrder',6, ...
         'HalfPowerFrequency1',100,'HalfPowerFrequency2',100+b, 'DesignMethod','butter', ...
         'SampleRate',14925);
bpFilt2 = designfilt('bandpassiir','FilterOrder',6, ...
         'HalfPowerFrequency1',100+b,'HalfPowerFrequency2',100+2*b, 'DesignMethod','butter', ...
         'SampleRate',14925);
fvtool(bpFilt, bpFilt2)
bpFilt2 = designfilt('bandpassfir','FilterOrder',1200, 'CutoffFrequency1',148,'CutoffFrequency2',152,...
    'SampleRate',14925);
fvtool(bpFilt2)

fs = 14925;
t = 0:1/fs:0.2;
x = [1]*sin(2*pi*[150 ]'.*t) ;
x = x(1:2048);
x2 = conv(bpFilt2.Coefficients, x); x2 = x2(1+1000:end-1000);
% x2 = filter(bpFilt, x);
figure;
X = fft(x); X2 = fft(x2);
subplot(2,1,1); plot(abs(X)); hold on; plot(abs(X2));
subplot(2,1,2); plot(angle(X)); hold on; plot(angle(X2));

figure; plot(x); hold on; plot(x2);
%%
recObj = audiorecorder(14762, 16, 1, -1);
recordblocking(recObj, 2);
y2 = getaudiodata(recObj);
% sound(y2, 14762);
spectrogram(y2,hamming(1024),512,4096, 14762, 'yaxis');