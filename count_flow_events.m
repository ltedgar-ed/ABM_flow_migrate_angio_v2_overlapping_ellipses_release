clear all
close all
clc

% curr_dir = pwd;
% save_directory = [pwd '\Simulation Data'];
% cd(save_directory)

root_file_name = 'ABM_output_ideal_cap_bed_3x3_Nt_192_krep_3.00_kmig_3.00_run';

runs = 25;

rev_events_runs = zeros(1,10);
drop_events_runs = zeros(1,10);

for r = 1:runs
    file_name = [root_file_name num2str(r) '.mat'];

    load(file_name)

    num_nodes = length(nodes);
    [num_vess num_timesteps] = size(vess_diameter);

    vess_conn = vess_conn + ones(num_vess, 2);

    % convert to uL/hr
    vess_flow = vess_flow/1e6;

    % convert to Pa
    nodal_pressures = nodal_pressures/12.96;

    rev_events = 0;
    drop_events = 0;
    restore_events = 0;

    flow_zero = 1e-5;

    for t = 2:num_timesteps
        for v = 1:num_vess

            if (vess_flow(v,t) < 0) && (vess_flow(v,t-1) > 0)
                rev_events = rev_events + 1;
            end

            if (vess_flow(v,t) > 0) && (vess_flow(v,t-1) < 0)
                rev_events = rev_events + 1;
            end

            if (abs(vess_flow(v,t)) < flow_zero) && (abs(vess_flow(v,t-1)) > flow_zero)
                drop_events = drop_events + 1;
            end

        end
    end

    rev_events_runs(r) = rev_events;
    drop_events_runs(r) = drop_events;
end

