close all
clc

%close(vidObj)

curr_dir = pwd;
save_directory = [pwd '\Simulation Data'];
cd(save_directory)
file_name = uigetfile;
load(file_name)

vidObj = VideoWriter([erase(file_name, '.mat') '_3D'],'MPEG-4');
vidObj.FrameRate = 30;
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

for t = 1:(plot_every/input.dt):num_timesteps
%for t = 1:1
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
        set(CYL,'FaceAlpha', 0.10)  
    end
    
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

        plot_ellipse_on_cylinder(Q, [x0; y0; 0], z_pos, c_pos, Z, R, Ain, Bin, alpha, 'r', true)

    end
    
    
    
    
    title([num2str(time(t)) ' hours '])
    
%     view(0, 20)
    
    axis off
    set(gca, 'FontSize', 24)
    set(gca, 'LineWidth', 2)
    set(figure(1), 'Color', 'w')
    
    fig = gcf;
    pos = fig.Position;
    set(fig, 'Position', [10 10 (pos(3)/pos(4))*700 700]);
    
    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)
end

close(vidObj);