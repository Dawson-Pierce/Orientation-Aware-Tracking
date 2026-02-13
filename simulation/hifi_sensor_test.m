%% High-Fidelity Data Generation
%
% Runs sim_quad for realistic quadrotor dynamics, resamples to sensor rate,
% and generates point cloud measurements.  Saves everything to
% hifi_data.mat for the tracking script (hifi_tracking.m).

clear; clc; close all

%% Import STL
TRsingle = stlread("drone_shape.stl");
TR = {TRsingle};

%% Define reference trajectory (circle -- easy to swap)
tf = 10;
dt_cmd = 0.01;
t_cmd = (0:dt_cmd:tf)';

r = 6;  w = 0.6;  center = [10, 15, 2];
p_cmd = [center(1) + r*cos(w*t_cmd), ...
         center(2) + r*sin(w*t_cmd), ...
         center(3) + (r*sin(w/2*t_cmd))];

%% Run high-fidelity sim
fprintf('Running sim_quad ... ');
[t_sim, p_sim, R_sim, ~, ~] = sim_quad(t_cmd, p_cmd);
fprintf('done  (%d steps)\n', numel(t_sim));

%% Plot comparison
figure;
plot3(p_cmd(:,1),p_cmd(:,2),p_cmd(:,3),'k'); hold on
plot3(p_sim(:,1),p_sim(:,2),p_sim(:,3),'b')
grid on; axis equal
xlabel('x'); ylabel('y'); zlabel('z')
legend('command','sim')
title('Trajectory comparison')

%% Resample to sensor rate
dt_sensor = 0.2;
t_sensor = (0:dt_sensor:tf)';
N = numel(t_sensor);

idx_sensor = zeros(N,1);
for k = 1:N
    [~, idx_sensor(k)] = min(abs(t_sim - t_sensor(k)));
end

p_sensor = p_sim(idx_sensor,:);
R_sensor = R_sim(:,:,idx_sensor);

%% Generate point cloud measurements
sensor_loc = [0 30 0];
sensor_dir = [0 -1 0.5]/sqrt(1.5);
sensor_slope = 0.01;

sen = sensor_station('location', sensor_loc, ...
    'view_vector', sensor_dir, 'grid_slope', sensor_slope);

points = cell(1, N);
fprintf('Generating point clouds ... ');
for k = 1:N
    Tk = {[R_sensor(:,:,k), p_sensor(k,:)'; 0 0 0 1]};
    points{k} = sen.get_points(TR, Tk);
end
fprintf('done\n');

%% Save
save('hifi_data.mat', ...
    't_sensor', 'dt_sensor', 'N', ...
    'p_sensor', 'R_sensor', 'points', ...
    'center', 'sensor_loc', 'sensor_dir', 'sensor_slope', ...
    'p_cmd', 'p_sim', 't_sim');

fprintf('Saved hifi_data.mat  (%d sensor frames)\n', N);


%% Plotting (optional)

plot_flag = 1;

if plot_flag
    figure; 
    h = scatter3(0,0,0,'.'); hold on
    traj = plot3(0,0,0,'k');

    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    for k = 1:length(points)
        temp_pts = points{k}; 
    
        h.XData = temp_pts(:,1);
        h.YData = temp_pts(:,2);
        h.ZData = temp_pts(:,3);

        traj.XData = p_sensor(1:k,1);
        traj.YData = p_sensor(1:k,2);
        traj.ZData = p_sensor(1:k,3);

        xlim([0 25])
        ylim([0 25])
        zlim([0 25])
    
        drawnow;
    
        pause(0.2)
    end
end