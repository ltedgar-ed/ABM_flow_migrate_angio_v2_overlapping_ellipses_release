close all
clear all
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

initial_cells = cells{1};

num_cells = length(initial_cells);

cell_vess_ID = zeros(num_cells, num_timesteps);

for t = 1:num_timesteps
    curr_cells = cells{t};
    
    [m n] = size(curr_cells);
    
    if (m ~= 0)
        for i = 1:m
            cell_ID = curr_cells(i,1);
            vess_ID = curr_cells(i,2);
        
            cell_vess_ID(cell_ID,t) = vess_ID;
        end
    end
end

