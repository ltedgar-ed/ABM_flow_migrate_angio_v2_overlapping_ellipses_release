close all
clc

%close(vidObj)

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

vidObj = VideoWriter([erase(file_name, '.mat') '_3D'],'MPEG-4');
vidObj.FrameRate = 20;
vidObj.Quality = 100;
open(vidObj);
plot_every = 0.5;

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

figure(1), hold on

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

load('ideal_cap_bed_3x3_flow_paths.mat')
[num_paths length_paths] = size(flow_paths_by_seg);

%flow_paths_by_seg = zeros(num_paths, length_paths-1);

% for p = 1:num_paths
%     path = flow_paths_by_nodes(p,:);
%     path = path + 1;
%     
%     for v = 1:length(path)-1
%         n0 = path(v);
%         n1 = path(v+1);
%         
%         node0 = nodes(n0,:);
%         node1 = nodes(n1,:);
%         
%         %plot([node0(1) node1(1)], [node0(2) node1(2)])
%         
%         for v2 = 1:num_vess
%             if (vess_conn(v2,1) == n0) && (vess_conn(v2,2) == n1)
%                 flow_paths_by_seg(p,v) = v2;
%             end
%         end
%     end
% end
% 
% for p = 1:num_paths
%     path = flow_paths_by_seg(p,:);
%     
%     for v = 1:length(path)
%         seg = path(v);
%         n0 = vess_conn(seg,1);
%         n1 = vess_conn(seg,2);
%         
%         node0 = nodes(n0,:);
%         node1 = nodes(n1,:);
%         
%         plot([node0(1) node1(1)], [node0(2) node1(2)], 'b')
%     end
% end

perfused_paths_over_time = [];
path_flow_over_time = [];

 for t = 1:num_timesteps
    perfused_paths = ones(num_paths,1);
    path_flow = zeros(num_paths,1);
    
    for p = 1:num_paths
        path = flow_paths_by_seg(p,:);

        for s = 1:length(path)
            seg = path(s);

            path_flow(p,1) = path_flow(p,1) + abs(vess_flow(seg,t));
            
            if (abs(vess_flow(seg,t)) < 1e-6) && (perfused_paths(p,1) == 1)
                perfused_paths(p,1) = 0;
            end
        end
    end

    perfused_paths_over_time = [perfused_paths_over_time perfused_paths];
    path_flow_over_time = [path_flow_over_time path_flow/sum(path_flow)];
    
    perfusion_loss_over_time(t) = 1 - sum(perfused_paths)/num_paths;
 end

perfusion_loss_over_time
%perfusion_loss_over_runs = [perfusion_loss_over_runs; perfusion_loss_over_time];
        
    