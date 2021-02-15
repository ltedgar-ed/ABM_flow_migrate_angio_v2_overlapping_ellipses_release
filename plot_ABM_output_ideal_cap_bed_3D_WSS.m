close all
clc

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

vidObj = VideoWriter([erase(file_name, '.mat') '_WSS_3D'],'MPEG-4');
vidObj.FrameRate = 30;
vidObj.Quality = 100;
open(vidObj);
plot_every = 1;

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

red_max = [0.2 0 0];
red_mid = [1 0 0];
red_min = [1.0 0.8 0.8];

WSS_max = max(max(vess_WSS));
WSS_max = 1.0;

shear_range = linspace(0, WSS_max, 8)';
leg_colormap = [];

for i = 1:length(shear_range)
    xi = shear_range(i)/WSS_max;
    
    if (xi < 0.5)
        leg_plot_color = 2*xi*red_mid + (1 - 2*xi)*red_min;
    else
        leg_plot_color = 2*(xi-0.5)*red_max + (1 - 2*(xi-0.5))*red_mid;
    end
    
    leg_colormap = [leg_colormap; leg_plot_color];
end

leg_tick_labels = cell(length(shear_range)+1,1);
tick_labels = linspace(0, WSS_max, length(shear_range)+1);

for i = 1:length(leg_tick_labels)
    leg_tick_labels{i} = num2str(round(tick_labels(i),1));
end

for t = 1:(plot_every/input.dt):num_timesteps
    curr_cells = cells{t};
    
    close (1)
    figure(1), hold on,
    axis([-0.1*maxmax 1.1*maxmax (ymax/2)-1.2*maxmax/2 (ymax/2)+1.2*maxmax/2 -1.2*maxmax/2 1.2*maxmax/2])
    
    for v = 1:num_vess
        x0 = nodes(vess_conn(v, 1), 1);
        y0 = nodes(vess_conn(v, 1), 2);
        x1 = nodes(vess_conn(v, 2), 1);
        y1 =  nodes(vess_conn(v, 2), 2);

        Z = norm([x1 - x0; y1 - y0]);
        R = vess_diameter(v,t)/2;
        
        xi = vess_WSS(v,t)/WSS_max;
    
        if (xi > 1.0)
            xi = 1.0;
        end
        
        if (xi < 0.0)
            xi = 0.0;
        end
        
        if (xi < 0.5)
            col = 2*xi*red_mid + (1 - 2*xi)*red_min;
        else
            col = 2*(xi-0.5)*red_max + (1 - 2*(xi-0.5))*red_mid;
        end
        
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
        
        if (vess_WSS(v,t) < 1e-3)
            CYL = surf(CXX, CYY, CZZ, 'FaceColor', [0.75 0.75 0.75]);
        else
            CYL = surf(CXX, CYY, CZZ, 'FaceColor', col);
        end
        
        set(CYL,'EdgeColor','None')
        set(CYL,'FaceAlpha', 0.75)  
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
    c.Label.String = ' WSS (Pa)';
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


