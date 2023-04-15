% Energy benchmarks

%% Energy revision
clear 
clc

% Edit these to select the correct data
n_steps = 10000;
benchmarks = ["Pool" "PoolCompSol" "Spring" "AngleBenchmark" "TorsionBenchmark" "Met" "T4LysozymeNoSolvent" "SolventBenchmark" "T4Lysozyme" ];
benchmark = benchmarks(7);
% ------------------------------------ %

workdir = "C:/PROJECTS/Quantom/Simulation/" + benchmark + "/Steps_" + string(n_steps)
file = fopen(workdir + "/energy.bin", "rb");
energy_data = fread(file, 'single');
fclose(file);


n_elements = length(energy_data)/3
energy_data = reshape(energy_data, [3, n_elements])';




%tiledlayout(2,1)
%nexttile

potE = energy_data(:,1);
kinE = energy_data(:,2);
totalE = energy_data(:,3);
%totalE = kinE + potE;
x = 1:length(potE);



plot(x, potE);
hold on
plot(x, kinE);
plot(x, totalE);
title(benchmark + " - Average energy")
legend("Potential energy", "Kinetic energy", "Total energy");
ylabel("Energy [J/mol]")
xlabel("time [fs]")
%xlim([0 100])
hold off
%label("Kinetic energy")


%
AN = 6.022*10^23;
m = 0.012 / AN
v = 1822;
k1 = 0.5 * 0.012 *v*v/AN
k2 = 0.5*m*v*v

