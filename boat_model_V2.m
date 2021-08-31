clear all; close all; clc;

s = tf('s');

%% Selected Parameters

% Boat dimensions
L = 11;                         % m, length of boat
W = 5;                          % m, width of boat
H = 2.5;                        % m, height of boat
rho = 640;                      % kg/m^3, density of boat (uniform)

% Initial conditions
x0 = [0, 0, 0];                 % initial position of [x, theta, phi]
v0 = [0, 0, 0];                 % initial d/dt of [x, theta, phi]

% System damping
R1 = 0.1;                       % Nms, Heave damping
R2 = 0.1;                       % Nms, Roll damping
R3 = 0.1;                       % Nms, Pitch damping

% Wave parameters
w = 1;                          % rad/s, wave frequency
max_height = .5;                % m, maximum height of wave
max_theta = 5*pi/180;           % rad, maximum theta angle of wave (roll)
max_phi = 10*pi/180;            % rad, maximum phi angle of wave (pitch)

% Simulation parameters
tf = 15;                        % s, time span of simulation
n = 900;                        % number of timesteps
lims = [10, 10, 10];            % x, y, and z axis limits (+/-)
fps = 30;                       % frame rate of animation video
filename = 'figures/Boat';      % name of generated mp4 file


%% Derived and Fixed Parameters

% Environmental constants
g = 9.81;       % m/s^2, gravity
rho_w = 1027;   % kg/m^3, saltwater density

% Boat parameters (https://dynref.engr.illinois.edu/rem.html)
m = W*L*H*rho;                  % kg, mass of boat
I_theta = (1/12)*m*(W^2+H^2);   % kg*m^2, roll moment of inertia
I_phi = (1/12)*m*(L^2+H^2);     % kg*m^2, pitch moment of inertia

% System spring constants
k1 = rho_w*g*W*L;             % kg/s, Heave spring constant
k2 = (1/12)*rho_w*g*L*W^3;    % Nm/rad, Roll spring constant
k3 = (1/12)*rho_w*g*W*L^3;    % Nm/rad, Pitch spring constant

% Vector forms
m_eq = [m, I_theta, I_phi];     % "mass-equivalent", i.e. - m or I
R = [R1, R2, R3];
k = [k1, k2, k3];
A = [max_height, max_theta, max_phi];

% Adjust initial conditions to make them absolute
x0 = x0 + [H/2-m*g/k1, -max_theta, -max_phi];

% Adjust limits vector for more convenient use
lims = [[-1, 1]; [-1, 1]; [-1, 1]] .* lims';

%% Laplace Domain

func = [w, s, s];

% Impedance matrix
Z1 = m*s + k1/s + R1;           % Heave impedance
Z2 = I_theta*s + k2/s + R2;     % Roll impedance
Z3 = I_phi*s + k3/s + R3;       % Pitch impedance
Z = [Z1, 0, 0;    
     0, Z2, 0;    
     0, 0, Z3 ];

% Effective forces
F0 = (k.*x0/s - m_eq.*v0)';         % Effects of initial conditions
F1 = max_height*k1*w/(s^2+w^2);     % Heave forcing
F1 = F1 - (rho_w*g*L*W*H/2-m*g)/s;  % Height equilibrium adjustment
F2 = max_theta*k2*s/(s^2+w^2);      % Roll forcing
F3 = max_phi*k3*s/(s^2+w^2);        % Pitch Forcing
F = F0 + [F1, F2, F3]';             % Comlete laplacian force vector

V = Z\F;
X = V/s + x0'/s;


%% Run Simulation

% Time frame
t = linspace(0, tf, n)';

% Generate response
h = sin(w*t)*A(1);          % Wave height h wrt time
theta = cos(w*t)*A(2);      % Wave angle theta wrt time
phi = cos(w*t)*A(3);        % Wave angle phi wrt time
wave = [h, theta, phi];
x = impulse(X, t);          % Position of boat wrt time

% Plot response
figure(1)
plot(t, x)
xlabel("Time [s]","FontSize",16)
ylabel("Position [m]/Angle [rad]","FontSize",16)
title("Absolute Boat Movement","FontSize",18)
legend(["Heave", "Roll", "Pitch"])
set(gcf, "PaperPosition", [0 0 5 5]);
set(gcf, "PaperSize", [5 5]);
saveas(gcf, "figures/Absolute.png");

figure(2)
plot(t, x-wave)
xlabel("Time [s]","FontSize",16)
ylabel("Position [m]/Angle [rad]","FontSize",16)
title("Boat Movement Relative to Wave","FontSize",18)
legend(["Heave", "Roll", "Pitch"])
set(gcf, "PaperPosition", [0 0 5 5]);
set(gcf, "PaperSize", [5 5]);
saveas(gcf, "figures/RelativeCombined.png");

names = ["Heave", "Roll", "Pitch"];
for indx=1:3
    figure(100+indx)
    plot(t, x(:,indx))
    title("Absolute Boat Movement: " + names(indx),"FontSize",18)
    xlabel("Time [s]","FontSize",16)
    if names(indx) == "Heave" 
        ylabel("Position [m]","FontSize",16)
    else
        ylabel("Angle [rad]","FontSize",16)
    end
    set(gcf, "PaperPosition", [0 0 5 5]);
    set(gcf, "PaperSize", [5 5]);
    saveas(gcf, "figures/Absolute" + names(indx) + ".png");
    
    figure(200+indx)
    plot(t, x(:,indx) - wave(:,indx))
    title("Boat Movement Relative to Wave: " + names(indx),"FontSize",18)
    xlabel("Time [s]","FontSize",16)
    if names(indx) == "Heave" 
        ylabel("Position [m]","FontSize",16)
    else
        ylabel("Angle [rad]","FontSize",16)
    end
    set(gcf, "PaperPosition", [0 0 5 5]);
    set(gcf, "PaperSize", [5 5]);
    saveas(gcf, "figures/Relative" + names(indx) + ".png");
end

%% Vibration Analysis
Y = inv(Z);                 % Admittance matrix
w = [0:0.001:1000];         % frequency vector

% Bode Plot
figure()
bode(Y/s, w);   % H = Y/s; Frequency response function

%% 2D Animation
frames = length(t);
x_max = max(x);

% Loop to create a 2D animation of each Heave, Roll, and Pitch
for indx=1:length(names)
    % Determine plot limits to fit everything nicely
    lim = max([x_max(indx), W/2, H/2, L/2]);
    lim = 1.1*lim;

    fig = figure("Name",names(indx) + "2D","NumberTitle","off");
    % https://www.mathworks.com/matlabcentral/answers/94495-how-can-i-create-animated-gif-images-in-matlab#answer_103847
    filename = "figures/" + names(indx) + '.gif';
    for frame=1:frames
        clf;
        if names(indx) == "Heave"
            draw_rectangle([0, x(frame,1)], W, H, 0);
            hold on
            draw_angled_line([0, wave(frame,1)], 3*lim, 0);
            hold off
        else
            draw_rectangle([0, 0], W, H, x(frame,indx));
            hold on
            draw_angled_line([0, 0], 3*lim, -wave(frame,indx));
            hold off
        end
        
        xlim([-lim,lim])
        ylim([-lim,lim])
        axis equal
        legend(["Frame: " + frame + "/" + frames])
        if names(indx) == "Heave"
            title(["Heave Animation " ...
                "COM Height: " + round(x(frame,1),2) + " m"])
        else
            title([names(indx) + " Animation " ...
                "Angle: " + round(x(frame,indx),2) + " rad"])
        end
        
        movieVector(frame) = getframe(fig, [0,0,560,420]);
    end

    myWriter = VideoWriter(filename, 'MPEG-4');
    myWriter.FrameRate = fps;

    open(myWriter);
    writeVideo(myWriter, movieVector);
    close(myWriter);
end

%% 3D Animation

fig = figure("Name","3D Animation","NumberTitle","off");

% Create animation slide for each timestep
for k=1:length(t)
    
    % Reset frame
    clf
    hold on
    
    % Plot surfaces of box (boat)
    [boat_vertices, boat_faces] = boatCoords(W, L, H, num2cell(x(k,:)));
    patch('Vertices', boat_vertices, 'Faces', boat_faces, 'FaceVertexCData',hsv(8),'FaceColor','flat')
    
    % Plot wave plane
    [X, Y, Z] = waveCoords(num2cell(wave(k,:)), lims);
    surf(X, Y, Z, 'FaceColor','b','FaceAlpha',.3,'EdgeAlpha',.3)
    
    % Format the plot
    grid on
    xlim(lims(1,:))
    ylim(lims(2,:))    
    zlim(lims(3,:))
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['t = ', num2str(t(k))])
    view([30 35])
    
    movieVector(k) = getframe(fig, [0,0,560,420]);
    
end

myWriter = VideoWriter(filename, 'MPEG-4');
myWriter.FrameRate = fps;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);


%% Animation Functions

% Return boat vertices and faces for use in patch function
function [vertices, faces] = boatCoords(W, L, H, p)

    % Boat position and orientation
    [h, theta, phi] = p{:};

    % Every n-length permutation of 1 and -1, where n=3
    % -- Never changes - can be moved outside of function/loop
    [a, b, c] = ndgrid([-1, 1]);
    direction = [a(:), b(:), c(:)];
    
    % Axis rotation matrices
    rotx = [1, 0, 0; 0, cos(phi), sin(phi); 0, -sin(phi), cos(phi)];
    roty = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
    
    % Vertex locations relative to boat orientation
    vertices = direction.*[L/2, W/2, H/2] + [0, 0, h];
    vertices = (roty*(rotx*vertices'))';
    
    % Identify faces by their vertices (based on directions list)
    % -- Never changes - can be moved outside of function/loop
    faces = [   1, 2, 4, 3;
                5, 6, 8, 7;
                1, 2, 6, 5;
                3, 4, 8, 7;
                1, 3, 7, 5;
                2, 4, 8, 6];
end

function [X, Y, Z] = waveCoords(q, lims)
    % Position and orientation of wave
    [h, theta, phi] = q{:};
    
    % Generate meshgrid plane
    xgrid = lims(1,1):1:lims(1,2);
    ygrid = lims(2,1):1:lims(2,2);
    [X,Y] = meshgrid(xgrid, ygrid);
    Z = zeros(size(X));
    plane(:,1,:) = X;
    plane(:,2,:) = Y;
    plane(:,3,:) = Z;
    
    % Axis rotation matrices
    rotx = [1, 0, 0; 0, cos(phi), sin(phi); 0, -sin(phi), cos(phi)];
    roty = [cos(theta), 0, sin(theta); 0, 1, 0; -sin(theta), 0, cos(theta)];
    
    % Rotate plane
    plane = pagemtimes(plane,rotx);
    plane = pagemtimes(plane,roty);
    
    % Return final X, Y and Z matrices
    X = squeeze(plane(:,1,:));
    Y = squeeze(plane(:,2,:));
    Z = squeeze(plane(:,3,:)) + h*ones(size(X));

end

% angled rectangle source
% https://www.mathworks.com/matlabcentral/answers/380257-how-do-i-draw-and-rotate-a-rectangle#comment_750053
function[]= draw_rectangle(center_location,width,height,theta)
    center1=center_location(1);
    center2=center_location(2);
    R = ([cos(theta), -sin(theta); sin(theta), cos(theta)]);
    X=([-width/2, width/2, width/2, -width/2]);
    Y=([-height/2, -height/2, height/2, height/2]);
    for i=1:4
        T(:,i)=R*[X(i); Y(i)];
    end
    x_lower_left=center1+T(1,1);
    x_lower_right=center1+T(1,2);
    x_upper_right=center1+T(1,3);
    x_upper_left=center1+T(1,4);
    y_lower_left=center2+T(2,1);
    y_lower_right=center2+T(2,2);
    y_upper_right=center2+T(2,3);
    y_upper_left=center2+T(2,4);
    x_coor=[x_lower_left x_lower_right x_upper_right x_upper_left];
    y_coor=[y_lower_left y_lower_right y_upper_right y_upper_left];
    patch('Vertices',[x_coor; y_coor]','Faces',[1 2 3 4],'Facecolor','none','Linewidth',1.2);
    axis equal;
end

function []=draw_angled_line(center_location, x_length, theta)
    x_center = center_location(1);
    y_center = center_location(2);
    x1 = x_center-x_length/2;
    x2 = x_center+x_length/2;
    y_length = x_length*tan(theta);
    y1 = y_center-y_length/2;
    y2 = y_center+y_length/2;
    plot([x1, x2], [y1, y2]);
end