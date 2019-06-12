function [b, newant_pos, new_ff] = Operate_mod(bb, ff, ant_pos, Num_ant, Num_fre, select_fre, select_ant, state_ant)

if nargin <= 7
    state_ant = 'same';
end
if nargin <= 6
    select_ant = 'uniform';
end
if nargin <= 5
    select_fre = 'uniform';
end

[N M] = size(bb);
b = zeros(Num_ant, Num_fre);

if strcmp(select_ant, 'uniform')
    ant_index  = ceil(linspace(1, N, Num_ant));
else
    ant_index  = ceil(N*rand(1, Num_ant));
end
ant_b = bb(ant_index, :);

newant_pos = ant_pos(ant_index, :);

if strcmp(select_fre, 'uniform')
    fre_index  = ceil(linspace(1, M, Num_fre));
    b = ant_b(:, fre_index);
    new_ff = ff(fre_index);
    new_ff = repmat(new_ff, Num_ant, 1);
elseif strcmp(select_fre, 'rand')
    if strcmp(state_ant, 'same')
        fre_index = ceil(M * rand(1, Num_fre));
        b = ant_b(:, fre_index);
        new_ff = ff(fre_index);
        new_ff = repmat(new_ff, Num_ant, 1);
    else
        fre_index = ceil(M * rand(Num_ant, Num_fre));
        b = ant_b(fre_index);
        new_ff = ff(fre_index);
    end
end

end
