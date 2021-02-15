close all
clc

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

vidObj = VideoWriter([erase(file_name, '.mat') '_3D'],'MPEG-4');
vidObj.FrameRate = 3;
vidObj.Quality = 100;
open(vidObj);

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

time = linspace(1,num_timesteps,num_timesteps) - 1;
time = time*t_time/24;

figure(1)

max_cell_num = 0;

for t = 1:num_timesteps
    if (max(cells{t}(:,1)) > max_cell_num)
        max_cell_num = max(cells{t}(:,1));
    end  
end

max_cell_speed_x = zeros(1,num_timesteps);
max_cell_speed_y = zeros(1,num_timesteps);

for t = 1:num_timesteps
    curr_cells = cells{t};
    
    for c = 1:length(curr_cells)
        cell_speed = curr_cells(c,9:10);
        
        if (-cell_speed(1) > max_cell_speed_x(t))
            max_cell_speed_x(t) = -cell_speed(1);
        end
        
        if (cell_speed(2) > max_cell_speed_y(t))
            max_cell_speed_y(t) = cell_speed(2);
        end
    end
end

red_max = [1 0.5 0.5];
red_min = [0.9 0 0];
blue_max = [0.5 0.5 1];
blue_min = [0 0 0.9];

cell_colors = zeros(max_cell_num,3);

for t = 1:num_timesteps
    curr_cells = cells{t};
    
    close (1)
    figure(1), hold on,
    axis([-50 150 -50 150 -100 100])
    
    for v = 1:num_vess
        x0 = nodes(vess_conn(v, 1), 1);
        y0 = nodes(vess_conn(v, 1), 2);
        x1 = nodes(vess_conn(v, 2), 1);
        y1 =  nodes(vess_conn(v, 2), 2);

        Z = norm([x1 - x0; y1 - y0]);
        R = vess_diameter(v,t)/2;
        
        [CYY CZZ CXX] = cylinder(R);
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
        set(CYL,'FaceAlpha', 0.25)  
    end
    
    for c = 1:length(curr_cells)
        id = curr_cells(c,1);
        vid = curr_cells(c,2);
        xi = curr_cells(c,3);
        zeta = curr_cells(c,4);
        pol = [curr_cells(c,5); curr_cells(c,6)];

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

        plot_ellipse_on_cylinder(Q, [x0; y0; 0], z_pos, c_pos, Z, R, Ain, Bin, alpha, 'r', true)

    end
    
    
    
    
    title([num2str(t-1) ' hours '])
    
    %view(-60, 30)
    axis off
    set(gca, 'FontSize', 24)
    set(gca, 'LineWidth', 2)
    set(figure(1), 'Color', 'w')
    
    fig = gcf;
    pos = fig.Position;
    set(fig, 'Position', [10 10 1.5*500 1.5*500]);

    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)
    
    pause(0.01)
end

close(vidObj);

curr_cells = cells{num_timesteps};

net_f_x = [];
net_f_y = [];

for c = 1:length(curr_cells)
    net_f_x = [net_f_x; curr_cells(c,7)];
    net_f_y = [net_f_y; curr_cells(c,8)];
end

vel_x = [];
vel_y = [];

for c = 1:length(curr_cells)
    vel_x = [vel_x; curr_cells(c,9)];
    vel_y = [vel_y; curr_cells(c,10)];
end


