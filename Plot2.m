clear all
close all


Files = dir('Data');
num_files = length(Files);
cd Data
data{:} = zeros(num_files-2);
for i = 3:num_files
   file = Files(i).name;
   data{i-2} = readmatrix(file);
end
cd ..
%% Video
% v = VideoWriter('Spikes.mp4', 'MPEG-4');
%v.LosslessCompression = true;
% v.FrameRate = 60;  % arbitrary
% open(v)
% f=figure;
% 
% pause(0.2) % let plot wake up
figure
hold on
plot_name = "Exponent Into Line: Viscosity = 0.05, cfl = 0.1";
xlim([0 2*pi])
xticks([0 pi 2*pi]);xticklabels({'0','\pi','2\pi'});
xlabel("X")
ylabel("Velocity")
title(plot_name)
for i=1:num_files - 2
    if  i == 1
        plot(data{i}(:,1), data{i}(:,2));
        data{i}(:,1);
        Legend_Name(1) = "Initial";
    elseif rem(i,200) == 0
        plot(data{i}(:,1), data{i}(:,2));
        data{i}(:,1);
        j = i/200 + 1;
        name = i*1 + " Iterations";
        Legend_Name(j) = name;
    end
    
end
legend(Legend_Name,'Location','best')
hold off