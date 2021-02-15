% plot_ellipse.m
% A function to add an ellipse to an existing 2D plot
% Lowell Taylor Edgar
% Usher Institute, University of Edinburgh
% 2019
% -------------------------------------------------------------------------
function [] = plot_ellipse_on_cylinder(ROT, x0, zcent, ccent, Z, R, A, B, Alpha, col, lw, fill_flag)

% Input arguments:
% xcent - x-coordinate of the ellipse centre
% ycent - y-coordinate of the ellipse centre
% A - radius along horizontal axis (prior to rotation)
% B - radius along vertical axis (prior to rotation)
% phi - rotation angle, counter-clockwise with respect of the horizontal axis (radians)
% col - color in which to print the ellipse

% A = Z*A;
% B = 2*pi*R*B;

% Number of subdivisions for plotting
sub_div = 100;

% Array of theta values (angle of rotation around centre)
phi = linspace(0, 2*pi, sub_div);

% Calculate the radius at each value of theta
for k = 1:sub_div
    r(k) = (A*B)/sqrt((B*cos(phi(k)))^2 + (A*sin(phi(k)))^2);
end

% Convert from polar to Cartersian coorindates
[z c] = pol2cart(phi, r);

line_div = sub_div;
zline = linspace(min(z), max(z), line_div);
cline = zeros(1,line_div);

zline2 = zeros(1,line_div);
cline2 = linspace(min(c), max(c), line_div);

zline3 = linspace(0, max(z), line_div/2);
cline3 = zeros(1,line_div/2);

% Create the rotation tensor
Q = [cos(Alpha) -sin(Alpha); sin(Alpha) cos(Alpha)];

% Rotate the ellipse by the angle alpha 
for k = 1:sub_div
    v = Q*[z(k); c(k)];
    z(k) = v(1);
    c(k) = v(2);
end

% Position the ellipse at the centre
z = z + zcent;
c = c + ccent;

% Rotate the ellipse by the angle alpha 
for k = 1:line_div
    v = Q*[zline(k); cline(k)];
    zline(k) = v(1);
    cline(k) = v(2);
end

zline = zline + zcent;
cline = cline + ccent;

% Rotate the ellipse by the angle alpha 
for k = 1:line_div
    v = Q*[zline2(k); cline2(k)];
    zline2(k) = v(1);
    cline2(k) = v(2);
end

zline2 = zline2 + zcent;
cline2 = cline2 + ccent;

% Rotate the ellipse by the angle alpha 
for k = 1:line_div/2
    v = Q*[zline3(k); cline3(k)];
    zline3(k) = v(1);
    cline3(k) = v(2);
end

zline3 = zline3 + zcent;
cline3 = cline3 + ccent;


[zcart, ycart, xcart] = pol2cart(c/R, R, z);
zcart = -zcart;
% [ycart2, zcart2, xcart2] = pol2cart(cline/R, R, zline);
% [ycart3, zcart3, xcart3] = pol2cart(cline2/R, R, zline2);

[zcart4, ycart4, xcart4] = pol2cart(cline3/R, R, zline3);
zcart4 = -zcart4;

for i = 1:length(xcart)
    vect_new = ROT*[xcart(i); ycart(i); zcart(i)] + x0;

    xcart(i) = vect_new(1);
    ycart(i) = vect_new(2);
    zcart(i) = vect_new(3);
end

for i = 1:length(xcart4)
    vect_new = ROT*[xcart4(i); ycart4(i); zcart4(i)] + x0;

    xcart4(i) = vect_new(1);
    ycart4(i) = vect_new(2);
    zcart4(i) = vect_new(3);
end



[zcartcent, ycartcent, xcartcent] = pol2cart(ccent/R, R, zcent);
zcartcent = -zcartcent;

vect_new = ROT*[xcartcent; ycartcent; zcartcent] + x0;

xcartcent = vect_new(1);
ycartcent = vect_new(2);
zcartcent = vect_new(3);

plot3(xcart, ycart, zcart, 'Color', col, 'LineWidth', lw), hold on,
% plot3(xcart2, ycart2, zcart2, 'w', 'LineWidth', 1.5)
% plot3(xcart3, ycart3, zcart3, 'w', 'LineWidth', 1.5)
plot3(xcart4, ycart4, zcart4, 'g', 'LineWidth', lw+0.25)

%plot3(xcartcent, ycartcent, zcartcent, 'r.', 'MarkerSize',30)