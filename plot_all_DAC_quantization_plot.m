clear;clc;
open('Figures/HW_impairments.fig');
h = gcf;

%% get data
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children'); 

objTypes = get(dataObjs{2,1}, 'Type');  %type of low-level graphics object


xdata_case3 = get(dataObjs{2,1}, 'XData');  %data from low-level grahics objects
ydata_case3 = get(dataObjs{2,1}, 'YData');

xdata_case2 = get(dataObjs{4,1}, 'XData');  %data from low-level grahics objects
ydata_case2 = get(dataObjs{4,1}, 'YData');

xdata_case1 = get(dataObjs{6,1}, 'XData');  %data from low-level grahics objects
ydata_case1 = get(dataObjs{6,1}, 'YData');
%%
figure(2);
hold on;
% Get the width and height of the figure
lbwh = get(1, 'position');
figw = lbwh(3);
figh = lbwh(4);
% Number of rows and columns of axes
ncols = 3;
nrows = 1;

% setup margin
marginw = 0.07;
marginh = 0.1;

% w and h of each axis in normalized units
axisw = (1 - (ncols+0.5) * marginw) / ncols;
axish = (1 / nrows) * 0.7
for ii=1:3
      % calculate the row and column of the subplot
      row = floor( ii/(ncols+1) ) + 1
      col = mod( ii-1, ncols ) + 1
      % calculate the left, bottom coordinate of this subplot
      axisl(ii) = (axisw+marginw) * (col-1) + marginw
      axisb(ii) = (axish+0.05) * (row-1) + 0.13
      %  plot the subplot
%       h= subplot('position', [axisl(ii), axisb(ii), axisw, axish] ); 
end


% color coding in figures
blue =  [0    0.4470    0.7410];
red = [0.8500    0.3250    0.0980];
yellow = [ 0.9290    0.6940    0.1250];
%%
figure(2)
h= subplot('position', [axisl(1), axisb(1), axisw, axish] );
set(0,'defaultTextInterpreter','latex'); %trying to set the default
% semilogx(DA(:,1),DA(:,2),'-o','linewidth',2);hold on
% semilogx(SA(:,1),SA(:,2),'-o','linewidth',2);hold on
% semilogx(FH(:,1),FH(:,2),'-o','linewidth',2);hold on
plot_order = [9,8,7,3,2,1,6,5,4];
gcf = plot(xdata_case1{plot_order(1),1},ydata_case1{plot_order(1),1},'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case1{plot_order(2),1},ydata_case1{plot_order(2),1},'--o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case1{plot_order(3),1},ydata_case1{plot_order(3),1},'-.o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case1{plot_order(4),1},ydata_case1{plot_order(4),1},'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(xdata_case1{plot_order(5),1},ydata_case1{plot_order(5),1},'--s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(xdata_case1{plot_order(6),1},ydata_case1{plot_order(6),1},'-.s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(xdata_case1{plot_order(7),1},ydata_case1{plot_order(7),1},'-x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(xdata_case1{plot_order(8),1},ydata_case1{plot_order(8),1},'--x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(xdata_case1{plot_order(9),1},ydata_case1{plot_order(9),1},'-.x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
grid on

% set(gca,'xtick',[64,128:128:1024])
set(gca,'FontSize', 14)

xlim([3,10])
ylim([40,62])
set(gca,'ytick',[40:5:62])
     
xlabel('DAC Quantization Number (Bits)')
ylabel('Spectral Efficiency (bps/Hz)')
% title('Dense Urban (LOS 100m, 58.8bps/Hz)')
t = title('Dense Urban MBB (Target SE 58.8bps/Hz)',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

legend('DA (U = 8)','DA (U = 16)','DA (U = 32)',...
       'SA (U = 8)', 'SA (U = 16)','SA (U = 32)',...
       'FH (U = 8)','FH (U = 16)','FH (U = 32)')
set(legend,'NumColumns',3);
set(legend,'Position',[axisl(1), axish+0.13, axisw, 0.15])

%%
h= subplot('position', [axisl(2), axisb(2), axisw, axish] );
set(0,'defaultTextInterpreter','latex'); %trying to set the default
% semilogx(DA(:,1),DA(:,2),'-o','linewidth',2);hold on
% semilogx(SA(:,1),SA(:,2),'-o','linewidth',2);hold on
% semilogx(FH(:,1),FH(:,2),'-o','linewidth',2);hold on
plot_order = [9,8,4,3,2,1,7,6,5];
gcf = plot(xdata_case2{plot_order(1),1},ydata_case2{plot_order(1),1},'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case2{plot_order(2),1},ydata_case2{plot_order(2),1},'--o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case2{plot_order(3),1},ydata_case2{plot_order(3),1},'-.o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case2{plot_order(4),1},ydata_case2{plot_order(4),1},'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(xdata_case2{plot_order(5),1},ydata_case2{plot_order(5),1},'--s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(xdata_case2{plot_order(6),1},ydata_case2{plot_order(6),1},'-.s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(xdata_case2{plot_order(7),1},ydata_case2{plot_order(7),1},'-x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(xdata_case2{plot_order(8),1},ydata_case2{plot_order(8),1},'--x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(xdata_case2{plot_order(9),1},ydata_case2{plot_order(9),1},'-.x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
grid on

% set(gca,'xtick',[64,128:128:1024])
set(gca,'FontSize', 14)

xlim([3,10])
ylim([3,5.2])
set(gca,'ytick',[3:0.5:5])
     
xlabel('DAC Quantization Number (Bits)')
ylabel('Spectral Efficiency (bps/Hz)')
% title('Dense Urban (LOS 100m, 58.8bps/Hz)')
t = title('50+Mbps Everywhere (Target SE 4.7bps/Hz)',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

legend('DA (U = 2)','DA (U = 4)','DA (U = 8)',...
       'SA (U = 2)', 'SA (U = 4)','SA (U = 8)',...
       'FH (U = 2)','FH (U = 4)','FH (U = 8)')
set(legend,'NumColumns',3);
set(legend,'Position',[axisl(2), axish+0.13, axisw, 0.15])
%%
h= subplot('position', [axisl(3), axisb(3), axisw, axish] );
set(0,'defaultTextInterpreter','latex'); %trying to set the default

plot_order = [2,1];
gcf = plot(xdata_case3{plot_order(1),1},ydata_case3{plot_order(1),1},'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(xdata_case3{plot_order(2),1},ydata_case3{plot_order(2),1},'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
grid on

% set(gca,'xtick',[64,128:128:1024])
set(gca,'FontSize', 14)

xlim([3,10])
ylim([3,13])
set(gca,'ytick',[3:2:13])
     
xlabel('DAC Quantization Number (Bits)')
ylabel('Spectral Efficiency (bps/Hz)')
% title('Dense Urban (LOS 100m, 58.8bps/Hz)')
t = title('Self-backhauling (Target SE 11.7bps/Hz)',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

legend('DA (U = 1)','SA and FH (U = 1)')
set(legend,'NumColumns',2);
set(legend,'Position',[axisl(3), axish+0.2, axisw, 0.05])

