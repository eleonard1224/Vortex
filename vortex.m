% Plots a vortex velocity field induced by a vortex tube in accordance
% with the Biot-Savart law

% Inputs
x = []; y = []; z = []; u = []; v = []; w = []; % x, y, z are
% spatial dimensions and u, v, w hold the velocity field
x1_x = 1; x1_y = 0; x1_z = 0; x2_x = 5; x2_y = 0; x2_z = 0; % These inputs
% hold the spatial coordinates for the two endpoints of the vortex tube
naxis = 20; nperp = 2; ncircular = 20; % These inputs hold the number of
% velocity vectors to plot along the axis of the vortex tube, perpendicular
% to the vortex tube, and circularly around the vortex tube
factor = -1/(4*pi);
sigma = 4;
gamma = 4e-04;
figurenumber = 1;



% Calculations
ihat = [1 0 0]; jhat = [0 1 0]; khat = [0 0 1]; % Unit vectors in the 
% x, y, z coordinate frame
x1 = [x1_x x1_y x1_z]; x2 = [x2_x x2_y x2_z]; % Coordinates of the two
% endpoints of the vortex tube
xhalf = 0.5 * (x1_x + x2_x); yhalf = 0.5 * (x1_y + x2_y); zhalf = 0.5 * (x1_z + x2_z);
shalf = [xhalf yhalf zhalf]; % Midpont of the vortex tube

% s, r, and t are axes which point along the vortex tube, 
% radially from the vortex tube, and tangentially around the vortex tube
% respectively
sx = x2_x - x1_x; sy = x2_y - x1_y; sz = x2_z - x1_z;
s = [sx sy sz]; smag = sqrt( sx*sx + sy*sy + sz*sz );
shat = s/smag;
if ( cross(shat, ihat) == 0 )
    r = cross(shat, jhat);
else
    r = cross(shat, ihat);
end
rx = r(1); ry = r(2); rz = r(3);
rmag = sqrt( rx*rx + ry*ry + rz*rz );
rhat = r/rmag;
t = cross(s, r);
tx = t(1); ty = t(2); tz = t(3);
tmag = sqrt( tx*tx + ty*ty + tz*tz );
that = t/tmag; 

% Calculating the velocity field - u, v, w 
for i = 1:naxis
    for j = 1:ncircular
        for k = 1:nperp
            sdist = (i-1)*smag/(naxis-1);
            theta = (j-1)*2*pi/ncircular;
            rdist = k * smag/(naxis-1);
            rpt = x1 + sdist*shat + rdist*cos(theta)*rhat + rdist*sin(theta)*that;
            x = [x rpt(1)]; y = [y rpt(2)]; z = [z rpt(3)];
            rvor = rpt - shalf;
            rvormag = sqrt( dot(rvor,rvor) );
            rvormagsigma = rvormag/sigma;
            rvormagsigma3 = rvormagsigma * rvormagsigma * rvormagsigma;
            phi = 1 - (1 - 1.5*rvormagsigma3)*exp(-rvormagsigma3); % phi approximates the local structure of the vortex tube
            vel = factor * gamma * phi * cross(rvor,s) / (rvormag^3);
            u = [u vel(1)]; v = [v vel(2)]; w = [w vel(3)];
        end
    end
end


% Plotting
figure(figurenumber)
quiver3(x, y, z, u, v, w, 'LineWidth', 2)
hold on
plot3([x1_x x2_x], [x1_y x2_y], [x1_z x2_z], 'g', 'LineWidth', 2)
set(gca,'FontSize',14)
grid off
set(gca,'visible','off')







