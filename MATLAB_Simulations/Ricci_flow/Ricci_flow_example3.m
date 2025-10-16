% Ricci_flow_example3
% Define the grid on the 2D plane
x = linspace(-10, 10, 100);
y = linspace(-10, 10, 100);
[X, Y] = meshgrid(x, y);

% Number of time steps and time step size
num_steps = 50;
time_step_size = 0.1;

% Initialize the metric tensor
metric = ones(size(X));

% Create a figure for animation
figure;
axis tight manual;
filename = 'ricci_flow_animation.gif';

% Evolve the metric over time steps and animate the process
for t = 1:num_steps
    % Ricci flow doesn't change the flat metric
    % But we can add a time-dependent factor for visualization purposes
    metric = metric + 0.1 * sin(t * time_step_size) * ones(size(X));
    
    % Plot the evolving 2D plane with updated metric
    surf(X, Y, metric);
    title(['2D Plane after Ricci Flow - Step ' num2str(t)]);
    xlabel('X');
    ylabel('Y');
    zlabel('Metric');
    
    % Adjust axis for better visualization
    axis([-10, 10, -10, 10, 0, 2]);
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
end
