clear all
close all


Files = dir('Data');
num_files = length(Files);
cd Data

for i = 3:num_files
      file = Files(i).name;
   data{i-2} = readmatrix(file);
end

https://www.scivision.dev/matlab-motion-jpeg-2000-mj2/#:~:text=Matlab%20can%20create%20lossless%20videos%20from%20image%20stacks,%28as%20described%20below%29%E2%80%93VideoWriter%20not%20available%20for%20GNU%20Octave.
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