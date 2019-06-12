function [I, xx, yy] = BPI(b, scene, imagesize, Ant, f0, df)
% function: back-propagation imaging with stepped frequency signal
% b: the measurement data, each row corresponds a receiving antenna
% scene = [Xmin, Xmax, Ymin, Ymax];the imaging domain
% Ant{1} the 3-D coordinates of the transmitting antennas
% Ant{2} the 3-D coordinates of the receiving antennas
% f0 the initial frequency
% df the frequency step

[M, N] = size(b);
c = 3e8; Ru = c / 2 / df;
Ni = 1024 * 32;
X = ifft(b.', Ni);
% imagesc(abs(X))
% construct the coordinate matrix of the image
Numx = imagesize(1); Numy = imagesize(2); I = zeros(Numy, Numx);
Xmin = scene(1); Xmax = scene(2); Ymin = scene(3); Ymax = scene(4);
Px = ones(Numy, 1)*linspace(Xmin, Xmax, Numx);
Py = linspace(Ymin, Ymax, Numy).'*ones(1, Numx);
Pz = ones(Numy, Numx) * 0;

Tx = Ant{1}(:, 1); Ty = Ant{1}(:, 2); Tz = Ant{1}(:, 3);
Rx = Ant{2}(:, 1); Ry = Ant{2}(:, 2); Rz = Ant{2}(:, 3);
%% %%%
Hwaitbar = waitbar(0, 'BP Imaging');
for kk = 1 : M
    R = sqrt((Tx(kk) - Px) .^ 2 + (Ty(kk) - Py) .^ 2 + (Tz(kk) - Pz) .^ 2)+...
        sqrt((Rx(kk) - Px) .^ 2 + (Ry(kk) - Py) .^ 2 + (Rz(kk) - Pz) .^ 2);
    r_prof = X(:, kk);
    R = mod(R, 2 * Ru);  	% if R > 2Ru, then it exceeds the unambiguous range
                      		% the aliasing effect is equivalent with mod(R, 2 * Ru)
    value = r_prof(ceil(R / (2 * Ru) * Ni));
    I = I + value .* exp(1j * 2 * pi * f0 * R / c);
    waitbar(kk / M, Hwaitbar);
end
close(Hwaitbar)
%%
xx = linspace(Xmin, Xmax, Numx);
yy = linspace(Ymin, Ymax, Numy);
end

