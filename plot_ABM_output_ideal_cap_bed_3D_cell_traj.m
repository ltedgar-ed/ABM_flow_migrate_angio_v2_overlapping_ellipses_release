clear all
close all
clc

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

% vidObj = VideoWriter([erase(file_name, '.mat') '_VELOCITY_3D'],'MPEG-4');
% vidObj.FrameRate = 15;
% vidObj.Quality = 100;
% open(vidObj);
% plot_every = 1.0;

cd(curr_dir)

%plot_cell_ID = 371;
%plot_cell_ID = 62;

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

% cell_vels = {};
% init_vels = zeros(length(cells{1}), 3);
% init_vels(:,1) = linspace(1, length(cells{1}), length(cells{1}));
% cell_vels{1} = init_vels;
no_near_neighbour_count = 0;

curr_cells = cells{1};
num_cells = length(curr_cells(:,1));

cell_traj_dist = zeros(num_cells,1);

for t = 2:num_timesteps
    curr_cells = cells{t};
    prev_cells = cells{t-1};
    
%     num_cells = length(curr_cells);
%     curr_vels = zeros(num_cells,3);
    
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
        
%         cell_vel = [d_Z d_C]/input.dt;
%         curr_vels(c,2:3) = cell_vel;
          cell_traj_dist(cell_ID,1) = cell_traj_dist(cell_ID,1) + abs(d_Z);
            
    end
end

[mean_rep_val mean_rep_ind] = min(abs(cell_traj_dist - mean(cell_traj_dist)));

plot_cell_ID1 = mean_rep_ind;

[min_rep_val min_rep_ind] = min(abs(cell_traj_dist));

plot_cell_ID2 = min_rep_ind;

[max_rep_val max_rep_ind] = max(abs(cell_traj_dist));
    
plot_cell_ID3 = max_rep_ind;

figure(1), hold on,

for v = 1:num_vess
    x0 = nodes(vess_conn(v, 1), 1);
    y0 = nodes(vess_conn(v, 1), 2);
    x1 = nodes(vess_conn(v, 2), 1);
    y1 =  nodes(vess_conn(v, 2), 2);

    Z = norm([x1 - x0; y1 - y0]);
    R = vess_diameter(v,1)/2;
    %R = 1;
    
    [CZZ CYY CXX] = cylinder(R);
    CZZ = -CZZ;
    CXX = Z*CXX;

    r = [x1 - x0; y1 - y0];
    x_vect = [1; 0];
    alpha = find_angle2D(x_vect, r);

    Q = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

    for j = 1:length(CXX)
        for i = 1:2
            xpt = CXX(i,j);
            ypt = CYY(i,j);
            zpt = CZZ(i,j);

            vect = [xpt; ypt; zpt];
            vect_new = Q*vect;

            CXX_new(i,j) = vect_new(1) + x0;
            CYY_new(i,j) = vect_new(2) + y0;
            CZZ_new(i,j) = vect_new(3);
        end
    end

    CXX = CXX_new;
    CYY = CYY_new;
    CZZ = CZZ_new;

    CYL = surf(CXX, CYY, CZZ); colormap(bone)
    set(CYL,'EdgeColor','None')
    set(CYL,'FaceAlpha', 0.05)  
end

    
for t = 1:num_timesteps
    curr_cells = cells{t};
    
%    close (1)
%    figure(1), hold on,
    axis([-0.1*maxmax 1.1*maxmax (ymax/2)-1.2*maxmax/2 (ymax/2)+1.2*maxmax/2 -1.2*maxmax/2 1.2*maxmax/2])    
    for c = 1:length(curr_cells)
        id = curr_cells(c,1);
        
        if (id == plot_cell_ID1) || (id == plot_cell_ID2) || (id == plot_cell_ID3)
            vid = curr_cells(c,2);
            xi = curr_cells(c,3);
            zeta = curr_cells(c,4);
            rho = curr_cells(c,5);
            pol = [curr_cells(c,6); curr_cells(c,7)];

            x0 = nodes(vess_conn(vid, 1), 1);
            y0 = nodes(vess_conn(vid, 1), 2);
            x1 = nodes(vess_conn(vid, 2), 1);
            y1 =  nodes(vess_conn(vid, 2), 2);

            Z = norm([x1 - x0; y1 - y0]);
            R = vess_diameter(vid,t)/2;

            r = [x1 - x0; y1 - y0];
            x_vect = [1; 0];
            alpha = find_angle2D(x_vect, r);

            Q = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];

            z_pos = Z*xi;
            c_pos = 2*pi*zeta*R;

            alpha = find_angle2D(pol, r);
            
            if (id == plot_cell_ID1)
                plot_color = [0 1 0];
                plot_ellipse_on_cylinder_color(Q, [x0; y0; 0], z_pos, c_pos, Z, R, Ain, Bin, alpha, plot_color, true)
            end
            
            if (id == plot_cell_ID2)
                plot_color = [1 0 0];
                plot_ellipse_on_cylinder_color(Q, [x0; y0; 0], z_pos, c_pos, Z, R, Ain, Bin, alpha, plot_color, true)
            end
            
            if (id == plot_cell_ID3)
                plot_color = [0 0 1];
                plot_ellipse_on_cylinder_color(Q, [x0; y0; 0], z_pos, c_pos, Z, R, Ain, Bin, alpha, plot_color, true)
            end
        end
    end
    
    %title([num2str(time(t)) ' hours '])
    
    axis off
    set(gca, 'FontSize', 24)
    set(gca, 'LineWidth', 2)
    set(figure(1), 'Color', 'w')
      
    fig = gcf;
    pos = fig.Position;
    set(fig, 'Position', [10 10 750 750]);
    
    % Write to the video file
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame)
end

mean_rep_ind
cell_traj_dist(mean_rep_ind)

min_rep_ind
cell_traj_dist(min_rep_ind)

max_rep_ind
cell_traj_dist(max_rep_ind)

%close(vidObj);