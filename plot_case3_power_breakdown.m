% This script is used to generate fig. 8 (block power breakdown) of the paper
% Exact data is from excel sheets.

clear;clc;

ASIC_DSP = 1;
DSP_scaling = 10;
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
DA_xdata_C3U1 = [128, 256, 384, 512];
DA_data_C3U1 = ...
[652.8	1305.6	1958.4	2611.2
106.4	97.9	92.9	89.4
2615.7	4448.9	6149.4	7785.2
1280.0	2560.0	3840.0	5120.0
7680.0	15360.0	23040.0	30720.0
0.0	0.0	0.0	0.0
5120.0	10240.0	15360.0	20480.0
164372.2	82381.2	54806.0	41098.7];

SA_xdata_C3U1 = [384,512,768,1024];
SA_data_C3U1 = ...
[0.0	0.0	0.0	0.0
168.9	168.9	168.9	168.9
6865.2	9153.1	13745.7	18346.1
480.0	640.0	960.0	1280.0
2880.0	3840.0	5760.0	7680.0
3840.0	5120.0	7680.0	10240.0
19200.0	25600.0	38400.0	51200.0
68364.1	51265.9	34263.2	25753.0];


FH_xdata_C3U1 = [384,512,768,1024];
FH_data_C3U1 = ...
[5.1	5.1	5.1	5.1
0.0	0.0	0.0	0.0
143.0	143.0	143.2	143.2
10.0	10.0	10.0	10.0
60.0	60.0	60.0	60.0
3840.0	5120.0	7680.0	10240.0
25600.0	34133.3	51200.0	68266.7
68364.1	51265.9	34263.2	25693.8];

if ~ASIC_DSP
    DA_data_C3U1(1,:) = DA_data_C3U1(1,:)*DSP_scaling;
    SA_data_C3U1(1,:) = SA_data_C3U1(1,:)*DSP_scaling;
    FH_data_C3U1(1,:) = FH_data_C3U1(1,:)*DSP_scaling;

end

%%

text_loc = 120;
font_size = 14;
figure
% ------------- Subplot for all Arrays --------------
b = bar(1:11,[DA_data_C3U1(:,2:4).';zeros(1,8);SA_data_C3U1(:,1:3).';zeros(1,8);FH_data_C3U1(:,1:3).']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,9,10,11])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',[' 256';' 384';' 512';' 384';' 512';' 768';' 384';' 512';' 768'])
xlim([0.5,11.5])
ylim([0,150])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')
% t = title('Digital Array',...
%     'Units', 'normalized',...
%     'Position', [0.5, 0.9, 0]); % Set Title with correct Position

text(2,text_loc,'DA','HorizontalAlignment','center','Fontsize',font_size)
text(6,text_loc,'SA','HorizontalAlignment','center','Fontsize',font_size)
text(10,text_loc,'FH','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Legends --------------
legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS',...
    'RF Amp','PA','Orientation','vertical')
set(legend,'Position',[0.4, axish+0.1, axisw, 0.05],'NumColumns',4)
