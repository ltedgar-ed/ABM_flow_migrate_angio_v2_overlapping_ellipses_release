% find_angle2D.m
% Find the angle between two 2D vectors and account for positive and
% negative angles
% Lowell Taylor Edgar
% Usher Institute, University of Edinburgh
% 2020
% -------------------------------------------------------------------------
function theta = find_angle2D(v1, v2)

v1 = v1/norm(v1);
v2 = v2/norm(v2);

v1dotv2 = dot(v1, v2);

if (v1dotv2 > 1.0)
    v1dotv2 = 1.0;
else
    if (v1dotv2 < -1.0)
        v1dotv2 = -1.0;
    end
end

theta = acos(v1dotv2);
theta_check = cross([v1; 0], [v2; 0]);

if (theta_check(3) < 0)
    theta = -theta;
end    