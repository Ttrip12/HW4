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
v = VideoWriter('test.mj2', 'Motion JPEG 2000');
v.LosslessCompression = true;
v.FrameRate = 5;  % arbitrary
open(v)
f=figure;

pause(0.2) % let plot wake up
for i=1:num_files - 2
    plot(data{i}(:,1), data{i}(:,2));
    plot_name = "Frame" + i;
    title(plot_name)
    xlim([0 2*pi])
    ylim([0 3.5])
    data{i}(:,1);
    im = getframe(f);
    writeVideo(v,im)
end
close(v)


% for i = 1:num_files
% figure
% hold on
% plot(data{i}(:,2))
% hold off
% end
% figure
% hold on
% ht = area(x,y,'facecolor',[.6 0 0]);
% for j = 1:num_files
%     x(2) = data{1}(1,1)
%     y(2) = data{1}(1,2)
%     set(ht,'XData',x)
%     set(ht,'YData',y)
%     drawnow
% end

% %
% 
% figure
% axes('XScale', 'linear', 'YScale', 'log')
% hold on
% title('10 Nodes, Second Order Accuracy Error, Second Order Derivative')
% loglog (der_2_2_err_10(4:13,1),abs(der_2_2_err_10(4:13,2)))
% xlabel('X Axis')
% ylabel('Y Axis')
% hold off
% 
% 
% 
% 
% 
% hold off
% function [A,B,C] = plotCSV3(file)
% 
% %     figure('Name', name);
% % 
% %     hold on
% 
%     A = readmatrix(append(file,'_.csv'));
%     B = readmatrix(append(file,'_100.csv'));
%     C = readmatrix(append(file,'_1000.csv'));
% 
% 
% 
% end
% 
% 
% function A = plotCSV(file, name)
% 
%     %figure('Name', name);
% 
%     A = readmatrix(file);
%     
%     %plot (A(:,1),A(:,2));
% 
% 
% end