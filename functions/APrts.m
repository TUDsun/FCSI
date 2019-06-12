function [fk, NN, E1, E2, mL, Mr, EH1, EH2, nL, nL1 ,dtn, Nr, newb, xx, yy] = ...
    APrts(f_par, scene, Imagesize, Ant, b, Num_fre, select_fre, state_ant)
% f_par = [f0 fmax df Num_fre]
% scene = [Xmin Xmax Ymin Ymax]
% Imagesize = [Numx Numy];
% Ant represents the positions of transmitter(Ant{1})  and receivers(Ant{2})
% b represents the received signal
% Num_fre represents the number of chosen frequency points
if nargin <= 7
    state_ant = 'same';
end
if nargin <= 6
    select_fre = 'uniform';
end

c = 3e8;
f0 = f_par(1); fmax = f_par(2); df = f_par(3);
k0 = f0 / df; kmax = fmax / df;
ff = f0 : df : fmax; N0 = length(ff);

Tx = Ant{1}(:, 1); Ty = Ant{1}(:, 2); Tz = Ant{1}(:, 3);
Rx = Ant{2}(:, 1); Ry = Ant{2}(:, 2); Rz = Ant{2}(:, 3);

Numx    = Imagesize(1); Numy = Imagesize(2); Num = Numx * Numy;
Xmin    = scene(1); Xmax = scene(2); Ymin = scene(3); Ymax = scene(4);
xx      = linspace(Xmin, Xmax, Numx);
Ptx     = repmat(xx, 1, Numy);
yy      = linspace(Ymin, Ymax, Numy);
Pty     = reshape((yy.' * ones(1, Numx)).', 1, Num);
Ptz     = zeros(1, Num);
NN      = size(Ant{1}, 1);
%-----------------------select the frequency and samples randomly-------------------------------
newb = zeros(NN, Num_fre);
if strcmp(select_fre, 'uniform')
    fk0  = ceil(linspace(1, N0, Num_fre));
    newb = b(:, fk0);
%     new_ff = ff(fk0);
%     new_ff = repmat(new_ff, NN, 1);
    fk0 = repmat(fk0, NN, 1);
elseif strcmp(select_fre, 'rand')
    if strcmp(state_ant, 'same')
        fk0 = ceil(N0 * rand(1, Num_fre));
        newb = b(:, fk0);
%         new_ff = ff(fk0);
%         new_ff = repmat(new_ff, NN, 1);
        fk0 = repmat(fk0, NN, 1);
    else
        fk0 = ceil(N0 * rand(NN, Num_fre));
        for ii = 1:NN
            newb(ii, :) = b(ii, fk0(ii, :));
%             new_ff(ii, :) = ff(fk0(ii, :));
        end
    end
end
fk = fk0.'+k0-1;
newb = reshape(newb.', NN * Num_fre, 1);
%----------------------------A * x -------------------------------------
M = 2 * kmax+2;
R = 1.5;
Mr = R * M;
i = 1;
while(i<Mr)
    i = i * 2;
end
Mr = i;
N = length(Ptx);
Msp = 3;
Nsp = 3;
tau = pi * Msp / (M * M * R * (R-0.5));
tn = zeros(NN, N).';
mL = zeros(NN, N).';
E1 = zeros(NN, 2 * Msp * N).';
E2 = zeros(NN, (kmax+1)).';
Hwaitbar=waitbar(0, '¹¹ÔìÏ¡Êè»ù¾ØÕó');
% A = zeros(NN * Num_fre, Num);
for nn = 1:NN
%     f = new_ff(nn, :).'; 
    R0 = sqrt((Tx(nn)-Ptx).^2+...
        (Ty(nn)-Pty).^2+...
        (Tz(nn)-Ptz).^2);
    R1 = sqrt((Rx(nn)-Ptx).^2+...
        (Ry(nn)-Pty).^2+...
        (Rz(nn)-Ptz).^2);
%     taul  = (R0+R1) / c;
%     index = (1:Num_fre)+(nn-1) * Num_fre;
%     A(index, :) = exp(-1j * 2 * pi * f * taul);
    
    tn(:, nn)  = 2 * pi * df * (R0+R1) / c;
    mL(:, nn) = floor(tn(:, nn) / 2 / pi * Mr);
    dt = tn(:, nn)-2 * pi * mL(:, nn) / Mr;
    tmp = exp(-(ones(2 * Msp, 1) * dt.'-2 * pi * (1-Msp:Msp).' * ones(1, N) / Mr).^2 / 4 / tau);
    E1(:, nn) = reshape(tmp, 2 * Msp * N, 1);
    E2(:, nn) = sqrt(pi / tau) * exp(tau * (0:kmax).^2).';
    waitbar(nn / NN);
end
close(Hwaitbar);
clear tmp
%-----------------------------A^H * x-----------------------------------
% Nmax = 2 * (ceil(max(tn(:) / 2 / pi * N))+1);
tau = pi * Nsp / (N * N * R * (R-0.5));
Nr = R * N;
i = 1;
while(i<Nr)
    i = i * 2;
end
Nr = i;
Ntn = ceil(max(tn(:) / 2 / pi * Nr))+1;
nL = zeros(NN, Num_fre).';
EH1 = zeros(NN, 2 * Nsp * Num_fre).';
EH2 = zeros(NN, Ntn).';
% fkp = fk / N * 2 * pi;
for nn = 1:NN
%     nL(:, nn) = floor(fkp(:, nn) / 2 / pi * Nr);
%     dt = fkp(:, nn)-2 * pi * nL(:, nn) / Nr;
%     tmp = exp(-(ones(2 * Nsp, 1) * dt.'-2 * pi * (1-Nsp:Nsp).' * ones(1, Num_fre) / Nr).^2 / 4 / tau);
    nL(:, nn) = fk(:, nn);
    tmp = exp(-(-2 * pi * (1-Nsp:Nsp).' * ones(1, Num_fre) / Nr).^2 / 4 / tau);
    EH1(:, nn) = reshape(tmp, 2 * Nsp * Num_fre, 1);
    EH2(:, nn) = sqrt(pi / tau) * exp(tau * (0:Ntn-1).^2).';
end
% nL1 = floor(tn / 2 / pi * Nr);
% dtn = tn / 2 / pi * Nr - nL1;
dtn = tn / 2 / pi * Nr;
nL1 = floor(dtn);%%%
dtn = dtn - nL1;%%%

end


