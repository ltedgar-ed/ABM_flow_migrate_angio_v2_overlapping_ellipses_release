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

time = linspace(0,num_timesteps,num_timesteps+1)*input.dt;

figure(1)

max_cell_num = 0;

for t = 1:num_timesteps
    if (max(cells{t}(:,1)) > max_cell_num)
        max_cell_num = max(cells{t}(:,1));
    end
end

rand_cell_pos = 2*(rand([max_cell_num,1])-0.5);

xmax = max(nodes(:,1));
ymax = max(nodes(:,2));

if (xmax > ymax)
    maxmax = xmax;
else
    maxmax = ymax;
end

node_degree = zeros(num_nodes,1);

for v = 1:num_vess
    node1 = vess_conn(v,1);
    node2 = vess_conn(v,2);
    node_degree(node1) = node_degree(node1) + 1;
    node_degree(node2) = node_degree(node2) + 1;
end

new_vess_conn = vess_conn;

right_edge_node1 = find(vess_conn == 155);
new_vess_conn(right_edge_node1) = 1;
right_edge_node2 = find(vess_conn == 174);
new_vess_conn(right_edge_node2) = 36;
right_edge_node3 = find(vess_conn == 198);
new_vess_conn(right_edge_node3) = 65;

plex_graph = graph;

for v = 1:num_vess
    plex_graph = addedge(plex_graph, new_vess_conn(v,1), new_vess_conn(v,2));
end

%plot(plex_graph)

cell_vels = {};
init_vels = zeros(length(cells{1}), 3);
init_vels(:,1) = linspace(1, length(cells{1}), length(cells{1}));
cell_vels{1} = init_vels;
no_near_neighbour_count = 0;

for t = 2:num_timesteps
    curr_cells = cells{t};
    prev_cells = cells{t-1};
    
    num_cells = length(curr_cells);
    curr_vels = zeros(num_cells,3);
              
    for c = 1:num_cells
        cell_ID = curr_cells(c,1);
        vess_ID = curr_cells(c,2);
        
        curr_vels(c,1) = cell_ID;
        x0 = nodes(vess_conn(vess_ID, 1), 1);
        y0 = nodes(vess_conn(vess_ID, 1), 2);
        x1 = nodes(vess_conn(vess_ID, 2), 1);
        y1 =  nodes(vess_conn(vess_ID, 2), 2);

        Z2 = norm([x1 - x0; y1 - y0]);
        R2 = vess_diameter(vess_ID,t)/2;
    
        index = find(prev_cells(:,1) == cell_ID);
        prev_vess_ID = prev_cells(index,2);
        
        d_Z = 0;
        d_C = 0;
        
        if (prev_vess_ID == vess_ID)
            d_xi = curr_cells(c,3) - prev_cells(index,3);
            d_zeta = curr_cells(c,4) - prev_cells(index,4);
            
            while (d_zeta > 1.0)
                d_zeta = d_zeta - 1.0;
            end
            
            while (d_zeta < 0.0)
                d_zeta = 1.0 + d_zeta;
            end
            
            d_Z = Z2*d_xi;
            d_C = R2*2*pi*d_zeta;
        end
        
        if (prev_vess_ID ~= vess_ID)
            xi1 = prev_cells(index,3);
            zeta1 = prev_cells(index,4);
            xi2 = curr_cells(c,3);
            zeta2 = curr_cells(c,4);
            
            x0 = nodes(vess_conn(prev_vess_ID, 1), 1);
            y0 = nodes(vess_conn(prev_vess_ID, 1), 2);
            x1 = nodes(vess_conn(prev_vess_ID, 2), 1);
            y1 =  nodes(vess_conn(prev_vess_ID, 2), 2);

            Z1 = norm([x1 - x0; y1 - y0]);
            R1 = vess_diameter(prev_vess_ID,t)/2;
            
            if (new_vess_conn(vess_ID,2) == new_vess_conn(prev_vess_ID,1))
                d_Z = -Z1*xi1 - Z2*(1-xi2);
            end
            
            if (new_vess_conn(vess_ID,1) == new_vess_conn(prev_vess_ID,2))
                d_Z = Z1*(1-xi1) + Z2*xi2;
            end
            
            if (new_vess_conn(vess_ID,1) == new_vess_conn(prev_vess_ID,1))
                d_Z = -Z1*xi1 + Z2*xi2;
            end
            
            if (new_vess_conn(vess_ID,2) == new_vess_conn(prev_vess_ID,2))
                d_Z = Z1*(1-xi1) - Z2*(1-xi2);
            end
            
            if (d_Z == 0)
                no_near_neighbour_count = no_near_neighbour_count + 1;
                d_Z = nan;
            end
             
            d_zeta = zeta2 - zeta1;

            while (d_zeta > 1.0)
                d_zeta = d_zeta - 1.0;
            end

            while (d_zeta < 0.0)
                d_zeta = 1.0 + d_zeta;
            end

            zeta1 = zeta2 - d_zeta;
            d_C = R2*2*pi*zeta2 - R1*2*pi*zeta1;

        end
        
        cell_vel = [d_Z d_C]/input.dt;
        curr_vels(c,2:3) = cell_vel;
    end
    
    cell_vels{t} = curr_vels;
end

no_near_neighbour_count

min_speed = 0;
max_speed = 0;

mean_speed = zeros(num_timesteps,1);
std_speed = zeros(num_timesteps,1);
all_speeds = [];

for t = 1:num_timesteps
    curr_cells = cells{t};
    curr_vels = cell_vels{t};
    mig_speeds = zeros(length(curr_vels),1);
    
    for c = 1:length(curr_vels)
        vess_ID = curr_cells(c,2);
        mig_vel_z = curr_vels(c,2);
        
        if (mig_vel_z > max_speed)
            max_speed = mig_vel_z;
        end    
        
        if (mig_vel_z < min_speed)
            min_speed = mig_vel_z;
        end
        
        mig_speeds(c) = sign(vess_flow(vess_ID,t))*mig_vel_z;
        
    end
    
    if (t ~= 1)
        all_speeds = [all_speeds; mig_speeds];
    end
    
    mean_speed(t) = mean(abs(mig_speeds));
    std_speed(t) = std(abs(mig_speeds));
end

figure(1), hist(all_speeds,2000)
xlabel(' mig speed (\mum/hr) ')
ylabel(' count ')
set(gcf,'Color','w')
set(gca,'FontSize',15)
set(gca,'LineWidth',1)
set(gca,'XLim',[-10 10])
set(gca,'YLim',[0 8000])