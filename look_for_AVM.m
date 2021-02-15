clear all
close all
clc

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

vidObj = VideoWriter([erase(file_name, '.mat') '_FLOW_3D'],'MPEG-4');
vidObj.FrameRate = 20;
vidObj.Quality = 100;
open(vidObj);

cd(curr_dir)

num_nodes = length(nodes);
[num_seg num_timesteps] = size(vess_diameter);
num_VESS = num_seg/5;

Pext = input.Pext;

vess_conn = vess_conn + ones(num_seg, 2);

% convert to uL/hr
vess_flow = vess_flow/1e6;

% convert to Pa
nodal_pressures = nodal_pressures/12.96;
transmural_pressures = nodal_pressures - Pext;

Ain = input.Ain;
Bin = input.Bin;

Lseg = 10;
cellspeed = 3;
t_time = Lseg/cellspeed;

time = linspace(0,num_timesteps,num_timesteps+1)*input.dt;

max_cell_num = 0;

for t = 1:num_timesteps
    if (max(cells{t}(:,1)) > max_cell_num)
        max_cell_num = max(cells{t}(:,1));
    end
end

mean_VESS_cells = zeros(num_VESS, num_timesteps);

for t = 1:num_timesteps
    for j = 0:num_VESS-1
        mean_VESS_cells(j+1,t) = mean(vess_num_cells(((j*5)+1):((j*5)+5),t));
    end
end

sum_VESS_cells = zeros(num_VESS, num_timesteps);

for t = 1:num_timesteps
    for j = 0:num_VESS-1
        sum_VESS_cells(j+1,t) = sum(vess_num_cells(((j*5)+1):((j*5)+5),t));
    end
end

VAR_sum_VESS_cells = var(sum_VESS_cells);

time_TEST = time(70:num_timesteps);
VAR_TEST = VAR_sum_VESS_cells(:, 70:num_timesteps);


