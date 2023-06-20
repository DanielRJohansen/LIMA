% Energy benchmarks

%% Energy revision
clear 
clc

% Edit these to select the correct data
n_steps = 3000;
benchmarks = ["Pool" "PoolCompSol" "Spring" "AngleBenchmark" "TorsionBenchmark" "Met" "T4LysozymeNoSolvent" "SolventBenchmark" "T4Lysozyme" "T4LysozymeNoSolventSmall"];
benchmark = "Spring";
%benchmark = benchmarks(10);
% ------------------------------------ %

workdir = "C:/PROJECTS/Quantom/Simulation/" + benchmark + "/Steps_" + string(n_steps)
file = fopen(workdir + "/energy.bin", "rb");
energy_data = fread(file, 'single');
fclose(file);

n_elements = length(energy_data)/3;
energy_data = reshape(energy_data, [3, n_elements])';

potE = energy_data(:,1);
kinE = energy_data(:,2);
totalE = energy_data(:,3);
%totalE = kinE + potE;

x = 1:length(potE);


from = 0;
to = inf;

% Plot original data
subplot(2,1,1)
plot(x, potE);
hold on
plot(x, kinE);
plot(x, totalE);
title(benchmark + " - Average energy")
legend("Potential energy", "Kinetic energy", "Total energy");
ylabel("Energy [J/mol]")
xlabel("time [fs]")
xlim([from to])
hold off

% Calculate and plot derivatives
dt = 1; % time step in fs
d_potE = diff(potE)/dt;
d_kinE = diff(kinE)/dt;
d_totalE = diff(totalE)/dt;
x_deriv = x(1:end-1);

subplot(2,1,2)
plot(x_deriv, d_potE);
hold on
plot(x_deriv, d_kinE);
plot(x_deriv, d_totalE);
title(benchmark + " - Energy derivatives")
legend("d_{pot}/dt", "d_{kin}/dt", "d_{total}/dt");
ylabel("Derivative of energy [J/mol/fs]")
xlabel("time [fs]")
xlim([from to])
hold off


%% 
%
AN = 6.022*10^23;
m = 0.012 / AN
v = 1822;
k1 = 0.5 * 0.012 *v*v/AN
k2 = 0.5*m*v*v

