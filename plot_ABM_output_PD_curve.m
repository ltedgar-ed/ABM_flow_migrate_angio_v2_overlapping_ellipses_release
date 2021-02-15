close all
clc

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

cd(curr_dir)

num_nodes = length(nodes);
[num_vess num_timesteps] = size(vess_diameter);

vess_conn = vess_conn + ones(num_vess, 2);

% convert to uL/hr
vess_flow = vess_flow/1e6;

% convert to Pa
nodal_pressures = nodal_pressures/12.96;

Ain = input.Ain;
Bin = input.Bin;

Lseg = 10;
cellspeed = 3;
t_time = Lseg/cellspeed;

time = linspace(0,num_timesteps,num_timesteps)*input.dt;

mean_diameter = mean(vess_diameter);

figure(1), plot(time, mean_diameter, 'LineWidth', 2.0)
xlabel(' time (hr) ')
ylabel(' diameter (\mum) ')

disp(['Equilbrium diameter of ' num2str(mean_diameter(num_timesteps))])
