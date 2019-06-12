%% Date: 19-April-2013
% Function: Back propagation imaging
% 'realdata_free.mat' is the real measurement data after system error correction 

clearvars; restoredefaultpath; 
addpath([pwd filesep 'spgl1-fast']);
addpath([pwd filesep 'functions']);

%% Select frequency components
load('realdata_free.mat');
findex = 1 : 2001;
b0 = b0(:, findex); b1 = b1(:, findex);
df = 1e6; ff = (findex + 499)*df; f0 = ff(1);
Num_ant  = 5; Num_fre  = 2001;

[b, ant_pos, new_ff] = Operate_mod(b1 - b0, ff, ant_pos, Num_ant, Num_fre, 'uniform', 'uniform', 'different');

%% BP Imaging
imagesize = [512, 512];         % set image pixels  
scene = [1, 5, 4, 8];           % set imaging domain              
% scene = [3.5, 4.5, 6, 7];              
TA = ant_pos; TA(:, 1) = ant_pos(:, 1) - 0.79/2; % 3-D coordinates of transmitting antennas
RA = ant_pos; RA(:, 1) = ant_pos(:, 1) + 0.79/2; % 3-D coordinates of receiving antennas
TA(:, 3) = 1.33; RA(:, 3) = 1.33;
Ant{1} = TA; Ant{2} = RA;       
[I, xx, yy] = BPI(b, scene, imagesize, Ant, f0, df);

%% show the imaging result
figure; 
imagesc(xx, yy, abs(I)); grid on; colorbar
title('BP Imaging', 'FontName',  'Times New Roman',  'FontSize',  10)
set(gca, 'YDir', 'normal',  'FontName',  'Times New Roman',  'FontSize',  10);
xlabel('Azimuth /m',  'FontName',  'Times New Roman',  'FontSize',  10)
ylabel('Range /m',  'FontName',  'Times New Roman',  'FontSize',  10)
axis equal tight
colormap jet
