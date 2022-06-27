scenarios = ["radar0.5m_source0.5m","radar0.5m_source1.0m","radar0.5m_source1.5m","radar0.5m_source2.0m"];
objects = ['layschip', 'smallcardboard', "smallpaperbag"];
lpf_cutoff = [1000, 800, 800];
% titles = ["paperbag", 'chips bag', 'cardboard'];
N=1; % questions per object scenario
qnum = 1;
for ii=1:length(objects)
    for jj=1:length(scenarios)
        permuted_digits = randperm(10)-1;
        for kk=1:N
%             answer = randi(10,1,9)-1;
            answer = permuted_digits((kk-1)*5+[1:9]);
%             if kk==1 % no repetition
%                 permuted_digits = randperm(10)-1;
%                 answer = permuted_digits(1:5);
%             else % allow repetition
%                 answer = randi(10,1,5)-1;
%             end           
            output_folder = sprintf('B:/listening_test/');
            output_filename = sprintf('%s/q%d.wav',output_folder, qnum);
            if ~exist(output_folder)
                mkdir(output_folder);
            end
            % concatenate wav files 
            y_out = [];
            for ll=1:9
                digit = answer(ll);
                wav_filename = sprintf('B:/experiment_data2/%s/%s/profile3/wav_eq/%d_%d.wav',...
                    objects(ii),scenarios(jj), digit,2);
                [y, fs] = audioread(wav_filename);
                 bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
                 'HalfPowerFrequency1',70,'HalfPowerFrequency2',lpf_cutoff(ii), 'DesignMethod','butter', ...
                 'SampleRate',fs);
                 y = filter(bpFilt, y);
                y_out = [y_out; 0.5*y(:)./max(abs(y)) ; 0.5*y(:)./max(abs(y)) ; 0.5*y(:)./max(abs(y)) ;];
            end
            audiowrite(output_filename, y_out, fs);
            save(sprintf('%s/q%d_%s_%s_%s.mat',output_folder, ...
                qnum, objects(ii), scenarios(jj),num2str(answer)), 'answer');
            qnum = qnum +1;
        end
    end
end

% y_out = [];
% for ll=0:9
%     wav_filename = sprintf('B:/experiment_data/microphone/wav_noeq_fs8k/%d_myvoice%d_%d.wav',...
%         ll,ll,2);
%     [y, fs] = audioread(wav_filename);
%     y_out = [y_out; 0.5*y(:)./max(abs(y)) ; 0.5*y(:)./max(abs(y)) ; ];
% end
% sound(0.1*y_out, fs);
% audiowrite('B:/listening_test/0-9.wav', y_out, fs);
%%
accuracy = [ 9 5 4 4 7;... %9
                      5 4 6 2 3;... %10
                      5 3 4 0 4;... %11
                      2 1 2 2 1;... %12
                      9 8 8 6 7;... %1
                      7 9 7 6 8;... %2
                      9 8 8 6 7;... %3
                      6 2 4 3 4;... %4
                      9 4 6 4 5;... %5
                      6 7 5 6 7;... %6
                      7 5 6 6 6;... %7
                      7 7 4 1 2;... %8
                      
                      ]./9; 
accuracy = accuracy(:, [1 2 3 5 ]);
mean_acc = reshape(mean(accuracy, 2), 4, []);
std_acc =  reshape(std(accuracy, 0,2), 4, []);
figure; 
for ii=1:size(mean_acc,2)
    e = errorbar(1:size(mean_acc,1), 100*mean_acc(:,ii), 100*std_acc(:,ii)/2); hold on;
    e.Marker = '*';
    e.MarkerSize = 10;
    e.LineWidth = 1.0;
%     e.Color = 'red';
%     e.CapSize = 15;
end
yticks(0:20:100);
xticklabels(arrayfun(@num2str, [0:20:100], 'UniformOutput', false));
ylabel('Recognition Accuracy (%)');
xticks(1:size(mean_acc,1));
xticklabels(arrayfun(@num2str, round(0.5*[1:size(mean_acc,1)],1), 'UniformOutput', false));
xlabel('Distanct (m)');
legend( ["paperbag", 'chips bag', 'cardboard']);

% Y = [4 5 2 3 7 1 0 8 9; 9 4 7 1 2 0 8 5 6;7 1 2 0 3 6 8 4 9;3 5 9 4 8 1 7 0 2];
% Y_hat = [];
% [c,order] = confusionmat(Y,Y_hat);
% cm = confusionchart(c); 
% cm.XLabel = 'Number';
% cm.Normalization = 'row-normalized';
%%
% filename = 'beatles_let_it_be_80db.wav';
filename = 'adele_someone_like_you_verse_80db.wav';
% filename = 'mccs0_sa1_80db.wav';
% filename = 'fadg0_sa1_80db.wav';
object = 'traderjoebag';
[y, fs] = audioread(sprintf('experiment_data2/%s/radar0.5m_source0.5m/profile3/wav_eq/%s', object, filename));
[y1, fs1] = audioread('audios/adele/someone_like_you_verse_60db.wav');
[llr, wss, stoi_val, nist_snr]=get_measures(y1, fs1, y, fs,0);
% [y, fs] = audioread(sprintf('experiment_data2/%s/radar0.5m_source1.0m/profile3/wav_eq/%s', object, filename));
% y = noiseReduction_YW(y, fs);
% y = y(fix(1*fs):fix(4.5*fs));
bpFilt = designfilt('bandpassiir','FilterOrder',20, ...
            'HalfPowerFrequency1',100,'HalfPowerFrequency2',1000, 'DesignMethod','butter', ...
            'SampleRate',fs);
y = filter(bpFilt, y);
yy = 0.05*y./max(abs(y));
sound(yy, fs);
clear sound 

audiowrite(sprintf('B:/demo/%s_%s',object,filename), yy, fs);
% fprintf("----------------------\n");
% for ii=0:9
%     for jj=1
%         %         [y, fs] = audioread(sprintf('experiment_data2/seaweed/radar0.5m_source0.5m/profile3/wav_eq/%d_%d.wav',ii,jj));
%         s = "radar0.5m_source2.0m";
%         [y, fs] = audioread(sprintf('B:/experiment_data2/smallpaperbag/%s/profile3/wav_eq/%d_%d.wav',s,ii,jj));
%         p1 = get_psnr(y,fs);
%         
% %         % denoise
% %         y = noiseReduction_YW(y, fs);
%         % filter
%         bpFilt = designfilt('bandpassiir','FilterOrder',10, ...
%             'HalfPowerFrequency1',130,'HalfPowerFrequency2',1500, 'DesignMethod','butter', ...
%             'SampleRate',fs);
%         y = filter(bpFilt, y);
%         
%         
%         p2 = get_psnr(y,fs);
%         fprintf("%.2f %.2f\n", p1,p2);
%         sound(0.3*y./max(abs(y)), fs);
%         pause(1.0);
%     end
% end