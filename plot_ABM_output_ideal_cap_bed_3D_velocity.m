close all
clc

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

vidObj = VideoWriter([erase(file_name, '.mat') '_velocity_3D'],'MPEG-4');
vidObj.FrameRate = 21;
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

cell_vels = {};
init_vels = zeros(length(cells{1}), 3);
init_vels(:,1) = linspace(1, length(cells{1}), length(cells{1}));
cell_vels{1} = init_vels;

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
        
        dZ = 0;
        dC = 0;
        
        if (prev_cells(index,2) == vess_ID)
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
        
        if (prev_cells(index,2) ~= vess_ID)
            prev_vess_ID = prev_cells(index,2);
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
            
            if (vess_conn(prev_vess_ID,2) == vess_ID)
                d_Z = Z1*(1 - xi1) + Z2*xi2;
            end
            
            if (vess_conn(prev_vess_ID,1) == vess_ID)
                d_Z = -Z1*xi1 - Z2*(1 - xi2) ;
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

min_speed = 0;
max_speed = 0;

mean_speed = zeros(num_timesteps,1);
std_speed = zeros(num_timesteps,1);

for t = 1:num_timesteps
    curr_vels = cell_vels{t};
    mig_speeds = zeros(length(curr_vels),1);
    
    for c = 1:length(curr_vels)
        mig_vel_z = curr_vels(c,2);
        
        if (mig_vel_z > max_speed)
            max_speed = mig_vel_z;
        end    
        
        if (mig_vel_z < min_speed)
            min_speed = mig_vel_z;
        end
        
        mig_speeds(c) = mig_vel_z;
    end
    
    mean_speed(t) = mean(abs(mig_speeds));
    std_speed(t) = std(abs(mig_speeds));
end

blue_max = [0 0.9 0];
blue_mid = [0.9 0.9 0.9];
blue_min = [0.9 0 0];

max_speed = 6;
min_speed = -6;

vel_width = max_speed - min_speed;
vel_range = linspace(min_speed, max_speed, 8)';
leg_colormap = [];

for i = 1:length(vel_range)
    xi = abs(vel_range(i) - min_speed)/vel_width;
    
    if (xi < 0.0)
        xi = 0.0;
    end
    
    if (xi > 1.0)
        xi = 1.0;
    end
        
    if (xi < 0.5)
        leg_plot_color = 2*xi*blue_mid + (1 - 2*xi)*blue_min;
    else
        leg_plot_color = 2*(xi-0.5)*blue_max + (1 - 2*(xi-0.5))*blue_mid;
    end
    
    leg_colormap = [leg_colormap; leg_plot_color];
end

leg_tick_labels = cell(length(vel_range)+1,1);
tick_labels = linspace(min_speed, max_speed, length(vel_range)+1);

for i = 1:length(leg_tick_labels)
    leg_tick_labels{i} = num2str(round(tick_labels(i),1));
end

for t = 1:(plot_every/input.dt):num_timesteps
    curr_cells = cells{t};
    curr_vels = cell_vels{t};
    
    close (1)
    figure(1), hold on,
    axis([-0.1*maxmax 1.1*maxmax (ymax/2)-1.2*maxmax/2 (ymax/2)+1.2*maxmax/2 -1.2*maxmax/2 1.2*maxmax/2])
    
%     for v = 1:num_vess
%         x0 = nodes(vess_conn(v, 1), 1);
%         y0 = nodes(vess_conn(v, 1), 2);
%         x1 = nodes(vess_conn(v, 2), 1);
%         y1 =  nodes(vess_conn(v, 2), 2);
% 
%         Z = norm([x1 - x0; y1 - y0]);
%         R = vess_diameter(v,t)/2;
%         
%         [CZZ CYY CXX] = cylinder(R);
%         CZZ = -CZZ;
%         CXX = Z*CXX;
%         
%         r = [x1 - x0; y1 - y0];
%         x_vect = [1; 0];
%         alpha = find_angle2D(x_vect, r);
%         
%         Q = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
% 
%         for j = 1:length(CXX)
%             for i = 1:2
%                 xpt = CXX(i,j);
%                 ypt = CYY(i,j);
%                 zpt = CZZ(i,j);
% 
%                 vect = [xpt; ypt; zpt];
%                 vect_new = Q*vect;
% 
%                 CXX_new(i,j) = vect_new(1) + x0;
%                 CYY_new(i,j) = vect_new(2) + y0;
%                 CZZ_new(i,j) = vect_new(3);
%             end
%         end
%         
%         CXX = CXX_new;
%         CYY = CYY_new;
%         CZZ = CZZ_new;
%         
%         CYL = surf(CXX, CYY, CZZ); colormap(bone)
%         set(CYL,'EdgeColor','None')
%         set(CYL,'FaceAlpha', 0.05)  
%     end
    
    for c = 1:length(curr_cells)
        id = curr_cells(c,1);
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

        index = find(curr_vels(:,1) == id);
        
        xi = abs(curr_vels(index,2) - min_speed)/vel_width;
    
        if (xi < 0.0)
            xi = 0.0;
        end

        if (xi > 1.0)
            xi = 1.0;
        end
    
        if (xi < 0.5)
            plot_color = 2*xi*blue_mid + (1 - 2*xi)*blue_min;
        else
            plot_color = 2*(xi-0.5)*blue_max + (1 - 2*(xi-0.5))*blue_mid;
        end
    
        plot_ellipse_on_cylinder_color(Q, [x0; y0; 0], z_pos, c_pos, Z, R, Ain, Bin, alpha, plot_color, true)
    end
    
    title([num2str(time(t)) ' hours '])
    
    axis off
    set(gca, 'FontSize', 24)
    set(gca, 'LineWidth', 2)
    set(figure(1), 'Color', 'w')
    
    c = colorbar;  
    colormap(leg_colormap);
    set(c, 'Ticks', linspace(c.Limits(1),c.Limits(2),9));
    set(c, 'TickLabels', leg_tick_labels);
    c.Label.String = ' migration velocity (\mum/hr)';
    c.Label.FontSize = 24;
    c.FontSize = 14;
    
    fig = gcf;
    pos = fig.Position;
    set(fig, 'Position', [10 10 750 750]);
    
    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)
end

close(vidObj);