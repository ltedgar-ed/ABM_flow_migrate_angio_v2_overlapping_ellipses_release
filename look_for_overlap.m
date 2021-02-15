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

time = linspace(0,num_timesteps,num_timesteps)*input.dt;

max_cell_num = 0;

for t = 1:num_timesteps
    if (max(cells{t}(:,1)) > max_cell_num)
        max_cell_num = max(cells{t}(:,1));
    end
end

mean_CELL_pos_overlap = zeros(1, num_timesteps);
mean_CELL_neg_overlap = zeros(1, num_timesteps);

for t = 1:num_timesteps
    curr_cells = cells{t};
    
    mean_CELL_pos_overlap(t) = mean(curr_cells(:,12));
    mean_CELL_neg_overlap(t) = mean(curr_cells(:,13));
end


figure(1), grid on, hold on, box on
plot(time, mean_CELL_pos_overlap, 'LineWidth', 1.0)
axis([0 100 0 1])
xlabel(' time (hours) ')
ylabel(' mean overlap among cells ')
set(gca, 'FontSize', 24)
set(gca, 'LineWidth', 2)
set(gcf, 'Color', 'w')
fig = gcf;
pos = fig.Position;
fig.Position = [1 2 1.5*pos(3) 1.5*pos(4)];

figure(2), grid on, hold on, box on
plot(time, mean_CELL_neg_overlap, 'LineWidth', 1.0)
axis([0 100 -20 0])
xlabel(' time (hours) ')
ylabel(' mean overlap among cells ')
set(gca, 'FontSize', 24)
set(gca, 'LineWidth', 2)
set(gcf, 'Color', 'w')
fig = gcf;
pos = fig.Position;
fig.Position = [1 2 1.5*pos(3) 1.5*pos(4)];
