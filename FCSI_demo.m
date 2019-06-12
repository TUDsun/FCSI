%% Date: 19-April-2013
% Function: process real measurement data with accelerated SPGL1
% 'realdata_free.mat' is the real measurement data after system error correction 

clearvars; restoredefaultpath; 
addpath([pwd filesep 'spgl1-fast']);
addpath([pwd filesep 'functions']);

    
%% select data
load('realdata_free.mat')
findex = 1 : 2001; df = 1e6; ff = (findex + 499) * df;
b0 = b0(:, findex); b1 = b1(:, findex);
f_par = [ff(1) ff(end) df];   % f_par = [f0 fmax df Num_fre]
imagesize = [512, 512];
scene = [1 5 4 8];
Num_fre = 250; select_fre = 'rand'; state_ant = 'different';
Num_ant = 5; 
aindex = ceil(linspace(1, size(b1, 1), Num_ant));
% aindex = ceil(rand(1, Num_ant)*size(b1, 1));
TA = ant_pos(aindex, :); TA(:, 1) = TA(:, 1) - 0.41;
RA = ant_pos(aindex, :); RA(:, 1) = RA(:, 1) + 0.41;
TA(:, 3) = 1.33; RA(:, 3) = 1.33;
Ant{1} = TA; Ant{2} = RA;
b = b1(aindex, :) - b0(aindex, :);
[fk, Na, E1, E2, mL, Mr, EH1, EH2, nL, nL1 ,dtn, Nr, bb, xx, yy] = ...
    APrts(f_par, scene, imagesize, Ant, b, Num_fre, select_fre, state_ant);
% [A,fk, Na, E1, E2, mL, Mr, EH1, EH2, nL, nL1 ,dtn, Nr, bb, xx, yy] = ...
%     APrts1(f_par, scene, imagesize, Ant, b, Num_fre, select_fre, state_ant);
% Na = the antenna number, fk -- Na * Num_fre; 
% E1 -- Na * 2 * Msp * N; the interpolation kernel
% E2 -- Na * (kmax + 1); 
% mL -- Na * N; the grid m that is closest to tn
% Mr = R * 2(kmax + 1); where, 2(kmax + 1) is the power of 2 

% EH1 -- Na * 2 * Nsp * Num_fre, the interpolation kernel  
% EH2 -- Na * Nmax / 2, 
% nL -- Na * Num_fre, represents the distance fk(:, nn)/N  
% nL1 -- Na * N
% dtn -- Na * N
% Nr = R * Nmax, where, Nmax is the power of 2 
%% SPG
sigma = 0.53 * norm(bb, 2);  % spg_bpdn needs to estimate the noise level             
options = spgSetParms('optTol', 3e-3, 'decTol', 5e-2);
% opts = spgSetParms('optTol', 1.5e-4); 
tau = 0; x = zeros(imagesize); x = x(:);
[x, r, g, info] = spgl1(fk, Na, E1, E2, mL, Mr, EH1, EH2, nL, nL1, dtn, Nr, bb, tau, sigma, x, options);
% [x,r,g,info] = spgl2(A, fk, Na, E1, E2, mL, Mr, EH1, EH2, nL, nL1, dtn, Nr, bb, tau, sigma, x, options);
x = reshape(x, imagesize);

%% plot image
figure
imagesc(xx, yy, abs(x.')); grid on; colorbar
ratio = roundn(Num_fre / (max(findex) - min(findex) + 1) * 100, -2); ratio = num2str(ratio);
set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 10);
title(['NUFFT SPG Imaging (',ratio, '% measurements)'],'FontName', 'Times New Roman', 'FontSize', 10)
xlabel('Azimuth /m', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('Range /m', 'FontName', 'Times New Roman', 'FontSize', 10)
colormap jet
axis equal tight

figure
plot(info.tauhis, info.rNorm2)
xlabel('\tau', 'FontName', 'Times New Roman', 'FontSize', 10)
ylabel('||r||_2', 'FontName', 'Times New Roman', 'FontSize', 10)
