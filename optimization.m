% lb = [-1,0,0,0,0];
% ub = [1,100,1000,1,2*pi];
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% amps = cand_amps;%[0.5486 3.491 0.8065 0.632]*1e6;
% dists = cand_dists; %[0.079 0.7032 0.8928 1.362];
% phis = ;
% lb = [zeros(1, length(amps)) max(d-0.044/2, 0) zeros(1, length(phis)) 0];
% ub = [Inf*ones(1,length(amps)) d+0.044/2 2*pi*ones(1, length(phis)) 2*pi];
% lb = [0 max(dists(2)-0.044, 0) 0];
% ub = [Inf dists(2)+0.044 2*pi];
% x0 = [amps d phis 1e-9];

for ii=1:8
    x0 = candidates(1:ii)
    lb = zeros(1, length(x0));
    ub = 2*pi*ones(1, length(x0));
    options = optimoptions('fmincon','Display','iter');
    % options=optimset('Algorithm','interior-point');
    tic;
    [x,fval,exitflag,output] = fmincon(@obj_func2,x0,[],[],Aeq,beq,lb,ub,[],options);
    toc;
    leftovernorm(ii) = fval;
    
end

fftsize = 512;
w = 2*pi*[0:fftsize-1].'/fftsize;
windowsize = 512;
A = exp(-1j.*(w-x).*(windowsize-1)/2).*diric((w-x), windowsize);
load('b.mat');
complex_gains = pinv(A)*b;
estimate = [x.' abs(complex_gains) angle(complex_gains)]
[candidates(1:8).' estimate(:,1)]
[estimate(:,2)] 
[estimate(:,3)]
ground_truth


%%

obj_func(x0)
obj_func(x)
mid = (length(x)-1)/3;
x(1:mid) - x0(1:mid)
x(mid+1:2*mid) - x0(mid+1:2*mid)
x(2*mid+1:3*mid) - x0(2*mid+1:3*mid)
% x0(1:4)
% x(1:4)
% x0(6:end)
% x(6:end)
% fprintf('AoA (deg): '); disp(rad2deg(acos(x(1:2:end))));
% fprintf('ToF (ns): '); disp(x(2:2:end));
% fprintf('Dop (Hz): '); disp(x(3:5:end));
% fprintf('Mag : '); disp(x(4:5:end));
% fprintf('Phase (rad): '); disp(x(5:5:end));

% clear cost;
% tofs = 0:100;
% aoa_deg = 1:180;
% dops = 0:100;
% for ii=1:length(dops)
%     for jj=1:length(tofs)
%         cost(ii, jj) = obj_func2([0, tofs(jj), dops(ii), 1, 0]);
%     end
% end
% figure;
% imagesc(cost);

%% test
f = [ 0.      ,  0.015625,  0.03125 ,  0.046875,  0.0625  ,  0.078125,...
        0.09375 ,  0.109375,  0.125   ,  0.140625,  0.15625 ,  0.171875,...
        0.1875  ,  0.203125,  0.21875 ,  0.234375,  0.25    ,  0.265625,...
        0.28125 ,  0.296875,  0.3125  ,  0.328125,  0.34375 ,  0.359375,...
        0.375   ,  0.390625,  0.40625 ,  0.421875,  0.4375  ,  0.453125,...
        0.46875 ,  0.484375, -0.5     , -0.484375, -0.46875 , -0.453125,...
       -0.4375  , -0.421875, -0.40625 , -0.390625, -0.375   , -0.359375,...
       -0.34375 , -0.328125, -0.3125  , -0.296875, -0.28125 , -0.265625,...
       -0.25    , -0.234375, -0.21875 , -0.203125, -0.1875  , -0.171875,...
       -0.15625 , -0.140625, -0.125   , -0.109375, -0.09375 , -0.078125,...
       -0.0625  , -0.046875, -0.03125 , -0.015625];
bw = 40e6;
f = f*bw + 5310e6;
lambs = 3e8 ./ f;
ant_spacing = lambs(1)/2;
tofs = [35e-9];
a = [0.8];
dop = [30];
aoa_cos = cos(deg2rad([ 78]));
N = 1; % rx_ant
H = zeros(64, N);
for i=1:length(tofs)
   tmp = a(i).* repmat(exp(-1j*2*pi*3e8*tofs(i)./lambs).',1,N) .* exp(-1j*2*pi*(ant_spacing*aoa_cos(i)./lambs.').*[0:N-1]);
   h = ifft(tmp, 64, 1);
   h = h .* repmat(exp(1j*2*pi*dop(i)/bw.*[0:63].'), 1, N);
   H = H + fft(h, 64, 1);
end
calibrated_H = H;
save calibrated_H.mat calibrated_H

% a= [0 1 0 0 0 0 0 0];
% H = fft(ifft(a .* exp(-1j*2*pi*d./lambs)) .* exp(1j*2*pi*dop/bw.*[0:63]), 64);
% rhs = a .* exp(-1j*2*pi*d./(3e8 ./ (f+dop)));
% figure(1);
% plot(angle(H)); hold on;
% plot(angle(rhs)); hold on;
