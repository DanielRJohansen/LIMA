%% Piecewise energies
clc
clear
workdir = "C:/PROJECTS/Quantom/Simulation/T4Lysozyme";
%workdir = "C:/PROJECTS/Quantom/Simulation/T4LysozymeNoSolventSmall";
file = fopen(workdir+"/PiecewiseEnergy.bin", "rb");
data = fread(file, 'single');
fclose(file);


% I want the potential energy to go now lower than 0, so i offset all the
% pots byt the largest negativ pot
minValue = min(data)
maxValue = max(data)
% Add the minimum value to every odd element
data(1:2:end) = data(1:2:end) + abs(minValue);

%max_step = 100;

% Reshape data into a matrix with n elements per column
%values_per_step = 192;
values_per_step = 25594;
data = reshape(data, values_per_step , []);

max_step = size(data, 2)

% Initialize the step variable
step = 1;

% Create a figure
figure;

% Create a slider
slider = uicontrol('Style', 'slider', 'Min', 1, 'Max', max_step, 'Value', step, ...
    'SliderStep', [1 / max_step, 10 / max_step], ...
    'Position', [100 50 200 20], 'Callback', {@sliderCallback, data, maxValue});

% Extract the first row of data
firstRow = data(:, step);

% Reshape the first row into pairs of elements
pairs = reshape(firstRow, 2, []);

% Create a stacked bar plot
bar(pairs', 'stacked');

% Slider callback function
function sliderCallback(source, ~, data, maxValue)
    % Get the slider value
    step = round(get(source, 'Value'));
    
    % Extract the corresponding row of data
    firstRow = data(:, step);
    
    % Reshape the row into pairs of elements
    pairs = reshape(firstRow, 2, []);
    
    % Update the stacked bar plot
    bar(pairs', 'stacked');
    %ylim([0, ceil(maxValue*2 / 1000) * 1000])
    ylim([0, 200000])
    %xlim([0 2000])
    title(sprintf('Step: %d', step));
end