% This script is used to generate fig. 8 (block power breakdown) of the paper
% Exact data is from excel sheets.

clear;clc;

ASIC_DSP = 1;
DSP_scaling = 3.2/(1/13);

grey = [0.83 0.82 0.78];

figure(1);
hold on;
% Get the width and height of the figure
lbwh = get(1, 'position');
figw = lbwh(3);
figh = lbwh(4);
% Number of rows and columns of axes
ncols = 3;
nrows = 1;

% setup margin
marginw = 0.05;
marginh = 0.1;

% w and h of each axis in normalized units
axisw = (1 - (ncols+0.5) * marginw) / ncols;
axish = (1 / nrows) * 0.75
for ii=1:3
      % calculate the row and column of the subplot
      row = floor( ii/(ncols+1) ) + 1
      col = mod( ii-1, ncols ) + 1
      % calculate the left, bottom coordinate of this subplot
      axisl(ii) = (axisw+marginw) * (col-1) + marginw
      axisb(ii) = (axish+0.05) * (row-1) + 0.125
      %  plot the subplot
%       h= subplot('position', [axisl(ii), axisb(ii), axisw, axish] ); 
end

%%
DA_xdata_U8 = 64:64:256;
DA_data_U8 = ...
[502.2	1004.3	1506.5	2008.6;
735.8	654.7	610.0	579.4;
1264.8	2132.5	2949.1	3739.5;
640.0	1280.0	1920.0	2560.0;
3840.0	7680.0	11520.0	15360.0;
0.0	0.0	0.0	0.0;
2560.0	5120.0	7680.0	10240.0;
68996.7	32198.0	20884.7	15446.4];

DA_xdata_U16 = 32,64:64:256;
DA_data_U16 = ...
[351.5	703.0	1406.0	2109.0	2812.1;
1523.8	1335.2	1164.3	1088.0	1088.0;
531.4	901.5	1618.3	2337.8	3117.1;
320.0	640.0	1280.0	1920.0	2560.0;
1920.0	3840.0	7680.0	11520.0	15360.0;
0.0	0.0	0.0	0.0	0.0;
1280.0	2560.0	5120.0	7680.0	10240.0;
31683.1	12157.1	5103.0	3212.4	2343.3];

DA_xdata_U32 = 32,64:64:256;
DA_data_U32 = ...
[552.4	1104.7	2209.5	3314.2	4419.0;
2791.0	2353.1	2176.0	2176.0	2176.0;
472.5	814.5	1558.5	2337.8	3117.1;
320.0	640.0	1280.0	1920.0	2560.0;
1920.0	3840.0	7680.0	11520.0	15360.0;
0.0	0.0	0.0	0.0	0.0;
1280.0	2560.0	5120.0	7680.0	10240.0;
16513.1	5430.4	2258.5	1402.3	1018.2];

SA_xdata = [256 512 768 1024];
SA_data_U8 = ...
[175.8	326.4	477.0	627.7
1091.4	1074.8	1066.8	1062.3
1453.9	2723.0	3960.3	5188.7
320.0	640.0	960.0	1280.0
1920.0	3840.0	5760.0	7680.0
2560.0	5120.0	7680.0	10240.0
12800.0	25600.0	38400.0	51200.0
225231.1	95066.7	58455.9	41862.8];

SA_data_U16 = ...
[251.1	401.7	552.4	703.0
1809.0	1702.4	1667.6	1652.3
757.2	1306.6	1874.9	2453.2
320.0	640.0	960.0	1280.0
1920.0	3840.0	5760.0	7680.0
2560.0	5120.0	7680.0	10240.0
12800.0	25600.0	38400.0	51200.0
134843.0	39249.0	21919.4	15210.0];

SA_xdata_U32 = [64,128,256,512];
SA_data_U32 = ...
[439.4	477.0	552.4	703.0
2176.0	2176.0	2351.3	2552.6
97.4	194.8	407.1	865.0
80.0	160.0	320.0	640.0
480.0	960.0	1920.0	3840.0
640.0	1280.0	2560.0	5120.0
3200.0	6400.0	12800.0	25600.0
28302.7	26596.7	21618.6	18023.0];


FH_xdata = 64:64:256;
FH_data_U8 = ...
[62.8	62.8	62.8	62.8
0.0	0.0	0.0	0.0
322.8	316.4	314.3	313.2
80.0	80.0	80.0	80.0
480.0	480.0	480.0	480.0
5120.0	10240.0	15360.0	20480.0
16213.3	32426.7	48640.0	64853.3
82571.1	39158.7	25634.7	19047.1];

FH_xdata_C1U16 = [32,64,128,192];
FH_data_C1U16 = ...
[175.8	175.8	175.8	175.8
0.0	0.0	0.0	0.0
342.1	309.8	294.8	290.7
160.0	160.0	160.0	160.0
960.0	960.0	960.0	960.0
5120.0	10240.0	20480.0	30720.0
14933.3	29866.7	59733.3	89600.0
46755.0	15842.7	6424.3	4025.6];

FH_xdata_C1U32 = [32,64,128,192];
FH_data_C1U32 = ...
[552.4	552.4	552.4	552.4
0.0	0.0	0.0	0.0
526.0	476.3	455.0	447.4
320.0	320.0	320.0	320.0
1920.0	1920.0	1920.0	1920.0
10240.0	20480.0	40960.0	61440.0
28586.7	57173.3	114346.7	171520.0
30048.9	8666.2	3234.7	1922.3];

if ~ASIC_DSP
    DA_data_U8(1,:) = DA_data_U8(1,:)*DSP_scaling;
    DA_data_U16(1,:) = DA_data_U16(1,:)*DSP_scaling;
    DA_data_U32(1,:) = DA_data_U32(1,:)*DSP_scaling;
    
    SA_data_U8(1,:) = SA_data_U8(1,:)*DSP_scaling;
    SA_data_U16(1,:) = SA_data_U16(1,:)*DSP_scaling;
    SA_data_U32(1,:) = SA_data_U32(1,:)*DSP_scaling;
    
    FH_data_U8(1,:) = FH_data_U8(1,:)*DSP_scaling;
    FH_data_C1U16(1,:) = FH_data_C1U16(1,:)*DSP_scaling;
    FH_data_C1U32(1,:) = FH_data_C1U32(1,:)*DSP_scaling;
end




%%
font_size = 14;
figure(1)

% ------------- Subplot for DA --------------
h= subplot('position', [axisl(1), axisb(1), axisw, axish] )
b = bar(1:11,[DA_data_U8(:,2:4).';zeros(1,8);DA_data_U16(:,1:3).';zeros(1,8);DA_data_U32(:,1:3).']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,9,10,11])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',['128';'192';'256';' 32';' 64';'128';...
                      ' 32';' 64';'128'])
xlim([0.5,11.5])
ylim([0,180])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')
t = title('Digital Array',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position
text(2,150,'U=8','HorizontalAlignment','center','Fontsize',font_size)
text(2,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(6,150,'U=16','HorizontalAlignment','center','Fontsize',font_size)
text(6,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(10,150,'U=32','HorizontalAlignment','center','Fontsize',font_size)
text(10,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
% legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS',...
%     'RF Amp','PA','Orientation','horizontal')
% set(legend,'Position',[0.4, axish+0.175, axisw, 0.05])
% ------------- Subplot for SA --------------

h= subplot('position', [axisl(2), axisb(2), axisw, axish] )
b = bar(1:11,[SA_data_U8(:,2:4).';zeros(1,8);SA_data_U16.';zeros(1,8);SA_data_U32(:,3:4).']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,8,10,11])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',[' 512';' 768';'1024';' 256';' 512';' 768';'1024';...
                      ' 256';' 512';])
xlim([0.5,11.5])
ylim([0,180])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')
t = title('Sub-Array',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position
text(2,150,'U=8','HorizontalAlignment','center','Fontsize',font_size)
text(2,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(6.5,150,'U=16','HorizontalAlignment','center','Fontsize',font_size)
text(6.5,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(10.5,150,'U=32','HorizontalAlignment','center','Fontsize',font_size)
text(10.5,140,'streams','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Subplot for FH --------------
h= subplot('position', [axisl(3), axisb(3), axisw, axish] )
b = bar(1:10,[FH_data_U8(:,1:3).';zeros(1,8);FH_data_C1U16(:,1:3).';zeros(1,8);FH_data_C1U32(:,1:2).']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,9,10])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',[' 64';'128';'192';' 32';' 64';'128';...
                      ' 32';' 64';'128'])
xlim([0.5,10.5])
ylim([0,180])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')
t = title('Fully-Connected Hybrid Array',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

text(2,150,'U=8','HorizontalAlignment','center','Fontsize',font_size)
text(2,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(6,150,'U=16','HorizontalAlignment','center','Fontsize',font_size)
text(6,140,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(9.5,150,'U=32','HorizontalAlignment','center','Fontsize',font_size)
text(9.5,140,'streams','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Legends --------------
legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS',...
    'RF Amp','PA','Orientation','horizontal')
set(legend,'Position',[0.4, axish+0.175, axisw, 0.05])
%% 4D Bar Figure Plot
% figure
% subplot(131)
% DA_data = zeros(3,4,8);
% DA_data(1,:,:) = DA_data_U8.';
% DA_data(2,:,:) = DA_data_U16.';
% DA_data(3,:,:) = DA_data_U32.';
% stacked_bar3(DA_data/1000);hold on
% title('Digital Array')
% colormap parula
% xticks([1,2,3,4])
% set(gca,'xticklabel',[' 64';'128';'192';'256'])
% yticks([1,2,3])
% set(gca,'yticklabel',[' 8';'16';'32'])
% grid on
% xlabel('Number of Array Elements')
% ylabel('Simultaneous Beams')
% zlabel('Power Consumption (W)')
% % legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS','RF Amp','PA')
% 
% subplot(132)
% SA_data = zeros(3,4,8);
% SA_data(1,:,:) = SA_data_U8.';
% SA_data(2,:,:) = SA_data_U16.';
% SA_data(3,:,:) = SA_data_U32.';
% stacked_bar3(SA_data/1000);hold on
% title('Sub-Array')
% colormap parula
% xticks([1,2,3,4])
% set(gca,'xticklabel',[' 256';' 512';' 768';'1024'])
% yticks([1,2,3])
% set(gca,'yticklabel',[' 8';'16';'32'])
% grid on
% xlabel('Number of Array Elements')
% ylabel('Simultaneous Beams')
% zlabel('Power Consumption (W)')
% % legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS','RF Amp','PA')
% 
% subplot(133)
% FH_data = zeros(3,4,8);
% FH_data(1,:,:) = FH_data_U8.';
% FH_data(2,:,:) = FH_data_U16.';
% FH_data(3,:,:) = FH_data_U32.';
% stacked_bar3(FH_data/1000);hold on
% title('Fully-Connected Hybrid Array')
% colormap parula
% xticks([1,2,3,4])
% set(gca,'xticklabel',[' 64';'128';'192';'256'])
% yticks([1,2,3])
% set(gca,'yticklabel',[' 8';'16';'32'])
% grid on
% xlabel('Number of Array Elements')
% ylabel('Simultaneous Beams')
% 
% zlabel('Power Consumption (W)')
% legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS',...
%     'RF Amp','PA','Orientation','horizontal')
