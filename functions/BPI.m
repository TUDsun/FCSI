function [I, xx, yy] = BPI(b,scene,imagesize,Ant,f0,df)
% 函数的功能是对步进频信号进行BP成像
% b表示接收数据，每个天线的接收数据按行排列
% scene = [Xmin, Xmax, Ymin, Ymax];成像场景的x，y坐标范围
% Ant{1}的每一行是发射天线的三维空间位置坐标，与数据b的接收顺序一致
% Ant{2}的每一行是接收天线的三维空间位置坐标，与数据b的接收顺序一致
% f0 表示起始频点
% df 表示频率步进量
[M,N] = size(b);
c = 3e8; Ru = c /2/df;
Ni = 1024*32;
X = ifft(b.',Ni);
% imagesc(abs(X))
% 构造图像位置坐标矩阵
Numx = imagesize(1); Numy = imagesize(2); I = zeros(Numy,Numx);
Xmin = scene(1); Xmax = scene(2); Ymin = scene(3); Ymax = scene(4);
Px = ones(Numy,1)*linspace(Xmin,Xmax,Numx);
Py = linspace(Ymin,Ymax,Numy).'*ones(1,Numx);
Pz = ones(Numy,Numx)*0;

Tx = Ant{1}(:,1); Ty = Ant{1}(:,2); Tz = Ant{1}(:,3);
Rx = Ant{2}(:,1); Ry = Ant{2}(:,2); Rz = Ant{2}(:,3);
%% %%%
Hwaitbar = waitbar(0,'BP 成像');
for kk = 1:M
    R = sqrt((Tx(kk)-Px).^2+(Ty(kk)-Py).^2+(Tz(kk)-Pz).^2)+...
        sqrt((Rx(kk)-Px).^2+(Ry(kk)-Py).^2+(Rz(kk)-Pz).^2);
    r_prof = X(:,kk);
    R = mod(R,2*Ru);  % 如果距离R>2Ru，说明成像方位超出了最大不模糊距离
                      % 场景的重叠相当于R对2Ru取余数运算——mod(R,2*Ru)
    value = r_prof(ceil(R/(2*Ru)*Ni));
    I = I+ value .* exp(1j*2*pi*f0*R/c);
    waitbar(kk/M,Hwaitbar);
end
close(Hwaitbar)
%%
xx = linspace(Xmin,Xmax,Numx);
yy = linspace(Ymin,Ymax,Numy);
end

