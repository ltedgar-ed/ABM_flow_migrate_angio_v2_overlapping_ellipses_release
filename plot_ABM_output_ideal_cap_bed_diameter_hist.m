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

all_diams = [];

for t = 1:num_timesteps
    all_diams = [all_diams; vess_diameter(:,t)];
end


figure(1), hist(all_diams,10)
xlabel(' vessel diameter (\mum) ')
ylabel(' count ')
set(gcf,'Color','w')
set(gca,'FontSize',15)
set(gca,'LineWidth',1)
set(gca,'XLim',[0 30])
% set(gca,'YLim',[0 8000])