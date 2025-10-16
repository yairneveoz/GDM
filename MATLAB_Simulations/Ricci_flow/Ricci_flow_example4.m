% Ricci_flow_example4
clear
clc
% Define the initial rectangular region
x = linspace(-5, 5, 50);
y = linspace(-3, 3, 30);
[X, Y] = meshgrid(x, y);

% Initialize the metric tensor with a non-uniform shape (initial distortion)
initial_metric = exp(-(X.^2 + Y.^2) / 5);

% Create a figure for animation
figure;
axis tight manual;
filename = 'ricci_flow_rectangle_animation.gif';

% Number of time steps and time step size
num_steps = 50;
time_step_size = 0.1;

% Evolve the metric over time steps and animate the process
for t = 1:num_steps
    % Apply Ricci flow to smooth out the metric
    smoothed_metric = initial_metric + 0.1 * del2(initial_metric) * time_step_size;
    
    % Plot the evolving rectangular region with updated metric
    surf(X, Y, smoothed_metric);
    title(['Rectangle after Ricci Flow - Step ' num2str(t)]);
    xlabel('X');
    ylabel('Y');
    zlabel('Metric');
    
    % Adjust axis for better visualization
    axis([-5, 5, -3, 3, 0, 1]);
    drawnow;
    
    % Capture the frame for animation
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);
    
    % Save the frame to the GIF file
    if t == 1
        imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime', 0.1);
    else
        imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
    end
    
    % Update initial metric for the next time step
    initial_metric = smoothed_metric;
end
