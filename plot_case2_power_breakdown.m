% This script is used to generate fig. 8 (block power breakdown) of the paper
% Exact data is from excel sheets.

clear;clc;

ASIC_DSP = 1;
DSP_scaling = 32;

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
DA_xdata_C2U2 = 384:128:1024;
DA_data_C2U2 = ...
[2109.0	2812.1	4218.1	5624.1
136.0	136.0	136.0	136.0
4675.6	6234.1	9351.2	12468.2
3840.0	5120.0	7680.0	10240.0
23040.0	30720.0	46080.0	61440.0
0.0	0.0	0.0	0.0
15360.0	20480.0	30720.0	40960.0
136718.8	102524.6	68206.9	51740.2];

DA_xdata_C2U4 = 256:128:768;
DA_data_C2U4 = ...
[1606.9	2410.3	3213.8	4820.7;
272.0	272.0	272.0	272.0;
3117.1	4675.6	6234.1	9351.2;
2560.0	3840.0	5120.0	7680.0;
15360.0	23040.0	30720.0	46080.0;
0.0	0.0	0.0	0.0;
10240.0	15360.0	20480.0	30720.0;
128182.4	85276.3	63801.1	42445.2];

DA_xdata_C2U8 = 256:128:768;
DA_data_C2U8 = ...
[2008.6	3012.9	4017.2	6025.8;
544.0	544.0	544.0	544.0;
3117.1	4675.6	6234.1	9351.2;
2560.0	3840.0	5120.0	7680.0;
15360.0	23040.0	30720.0	46080.0;
0.0	0.0	0.0	0.0;
10240.0	15360.0	20480.0	30720.0;
97236.3	65437.7	49411.5	33176.3];


SA_xdata_C2U2 = [768 1024 1536 2048];
SA_data_C2U2 = ...
[453.5	604.2	905.4	1206.7;
171.8	172.2	172.1	172.1;
1393.2	1862.2	2791.4	3721.8;
960.0	1280.0	1920.0	2560.0;
5760.0	7680.0	11520.0	15360.0;
7680.0	10240.0	15360.0	20480.0;
38400.0	51200.0	76800.0	102400.0;
168588.6	128477.9	85276.3	63948.2];

SA_xdata_C2U2 = [768 1024 1536 2048];
SA_data_C2U4 = ...
[458.2	608.9	910.2	1211.4;
285.3	286.3	286.7	286.6;
1199.2	1602.2	2405.1	3206.7;
960.0	1280.0	1920.0	2560.0;
5760.0	7680.0	11520.0	15360.0;
7680.0	10240.0	15360.0	20480.0;
38400.0	51200.0	76800.0	102400.0;
206456.4	158062.3	106127.6	79584.5];

SA_data_C2U8 = ...
[477.0	627.7	929.0	1230.3;
544.0	544.0	544.0	544.0;
1168.9	1558.5	2337.8	3117.1;
960.0	1280.0	1920.0	2560.0;
5760.0	7680.0	11520.0	15360.0;
7680.0	10240.0	15360.0	20480.0;
38400.0	51200.0	76800.0	102400.0;
299108.2	232181.9	162865.2	124402.3];

FH_xdata_C2U2 = [384,512,768,1024];
FH_data_C2U2 = ...
[11.0	11.0	11.0	11.0;
146.3	145.9	145.5	145.1;
60.6	60.6	60.7	60.7;
20.0	20.0	20.0	20.0;
120.0	120.0	120.0	120.0;
7680.0	10240.0	15360.0	20480.0;
35840.0	47786.7	71680.0	95573.3;
170148.6	127300.0	85354.9	64095.6];

FH_xdata_C2U4 = [384,512,768,1024];
FH_data_C2U4 = ...
[25.1	25.1	25.1	25.1;
235.5	235.0	234.7	234.3;
48.6	48.6	48.5	48.5;
40.0	40.0	40.0	40.0;
240.0	240.0	240.0	240.0;
15360.0	20480.0	30720.0	40960.0;
56320.0	75093.3	112640.0	150186.7;
106470.2	79767.9	52945.4	39612.1];

FH_xdata_C2U8 = [192,256,384,512];
FH_data_C2U8 = ...
[62.8	62.8	62.8	62.8
0.0	0.0	0.0	0.0
97.4	97.4	97.4	97.4
80.0	80.0	80.0	80.0
480.0	480.0	480.0	480.0
15360.0	20480.0	30720.0	40960.0
48640.0	64853.3	97280.0	129706.7
168588.6	126715.1	84884.5	63801.1];

if ~ASIC_DSP
    DA_data_C2U2(1,:) = DA_data_C2U2(1,:)*DSP_scaling;
    DA_data_C2U4(1,:) = DA_data_C2U4(1,:)*DSP_scaling;
    DA_data_C2U8(1,:) = DA_data_C2U8(1,:)*DSP_scaling;
    
    SA_data_C2U2(1,:) = SA_data_C2U2(1,:)*DSP_scaling;
    SA_data_C2U4(1,:) = SA_data_C2U4(1,:)*DSP_scaling;
    SA_data_C2U8(1,:) = SA_data_C2U8(1,:)*DSP_scaling;
    
    FH_data_C2U2(1,:) = FH_data_C2U2(1,:)*DSP_scaling;
    FH_data_C2U4(1,:) = FH_data_C2U4(1,:)*DSP_scaling;
    FH_data_C2U8(1,:) = FH_data_C2U8(1,:)*DSP_scaling;
end

%%

text_loc = 300;
font_size = 14;
figure
% ------------- Subplot for DA --------------
h= subplot('position', [axisl(1), axisb(1), axisw, axish] )
b = bar(1:12,[DA_data_C2U2(:,2:4).';zeros(1,8);DA_data_C2U4(:,2:4).';zeros(1,8);DA_data_C2U8.']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,9,10,11,12])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',[' 512';' 768';'1024';' 384';' 512';' 768';...
                      ' 256';' 384';' 512';' 768'])
xlim([0.5,12.5])
ylim([0,350])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')
t = title('Digital Array',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

text(2,text_loc,'U=2','HorizontalAlignment','center','Fontsize',font_size)
text(2,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(6.5,text_loc,'U=4','HorizontalAlignment','center','Fontsize',font_size)
text(6.5,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(10.5,text_loc,'U=8','HorizontalAlignment','center','Fontsize',font_size)
text(10.5,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Subplot for SA --------------

h= subplot('position', [axisl(2), axisb(2), axisw, axish] )
b = bar(1:10,[SA_data_C2U2(:,2:4).';zeros(1,8);SA_data_C2U4(:,2:4).';zeros(1,8);SA_data_C2U8(:,3:4).']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,9,10])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',['1024';'1536';'2048';'1024';'1536';'2048';...
                      '1536';'2048'])
xlim([0.5,10.5])
ylim([0,350])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')

t = title('Sub-Array',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

text(2,text_loc,'U=2','HorizontalAlignment','center','Fontsize',font_size)
text(2,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(6,text_loc,'U=4','HorizontalAlignment','center','Fontsize',font_size)
text(6,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(9.5,text_loc,'U=8','HorizontalAlignment','center','Fontsize',font_size)
text(9.5,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Subplot for FH --------------
h= subplot('position', [axisl(3), axisb(3), axisw, axish] )
b = bar(1:11,[FH_data_C2U2(:,2:4).';zeros(1,8);FH_data_C2U4(:,1:3).';zeros(1,8);FH_data_C2U8(:,1:3).']/1000,'stacked');
colormap parula
b(8).FaceColor = grey;
xticks([1,2,3,5,6,7,9,10,11])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',[' 512';' 768';'1024';' 384';' 512';' 768';...
                      ' 192';' 256';' 384'])
xlim([0.5,11.5])
ylim([0,350])
grid on
xlabel('Number of Array Elements')
ylabel('Power Consumption (W)')
t = title('Fully-Connected Hybrid Array',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

text(2,text_loc,'U=2','HorizontalAlignment','center','Fontsize',font_size)
text(2,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(6,text_loc,'U=4','HorizontalAlignment','center','Fontsize',font_size)
text(6,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)
text(10,text_loc,'U=8','HorizontalAlignment','center','Fontsize',font_size)
text(10,text_loc-20,'streams','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Legends --------------
legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS',...
    'RF Amp','PA','Orientation','horizontal')
set(legend,'Position',[0.4, axish+0.175, axisw, 0.05])
%%
% figure
% subplot(131)
% DA_data = zeros(3,4,8);
% DA_data(1,:,:) = DA_data_U2.';
% DA_data(2,:,:) = DA_data_U4.';
% DA_data(3,:,:) = DA_data_U8.';
% stacked_bar3(DA_data/1000);hold on
% title('Digital Array')
% colormap parula
% xticks([1,2,3,4])
% set(gca,'xticklabel',['256';'384';'512';'768'])
% yticks([1,2,3])
% set(gca,'yticklabel',['2';'4';'8'])
% grid on
% xlabel('Number of Array Elements')
% ylabel('Simultaneous Beams')
% zlabel('Power Consumption (W)')
% % legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS','RF Amp','PA')
% 
% subplot(132)
% SA_data = zeros(3,4,8);
% SA_data(1,:,:) = SA_data_U2.';
% SA_data(2,:,:) = SA_data_U4.';
% SA_data(3,:,:) = SA_data_U8.';
% stacked_bar3(SA_data/1000);hold on
% title('Sub-Array')
% colormap parula
% xticks([1,2,3,4])
% set(gca,'xticklabel',[' 256';' 512';' 768';'1024'])
% yticks([1,2,3])
% set(gca,'yticklabel',['2';'4';'8'])
% grid on
% xlabel('Number of Array Elements')
% ylabel('Simultaneous Beams')
% zlabel('Power Consumption (W)')
% % legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS','RF Amp','PA')
% 
% subplot(133)
% FH_data = zeros(3,4,8);
% FH_data(1,:,:) = FH_data_U2.';
% FH_data(2,:,:) = FH_data_U4.';
% FH_data(3,:,:) = FH_data_U8.';
% stacked_bar3(FH_data/1000);hold on
% title('Fully-Connected Hybrid Array')
% colormap parula
% xticks([1,2,3,4])
% set(gca,'xticklabel',['256';'384';'512';'768'])
% yticks([1,2,3])
% set(gca,'yticklabel',['2';'4';'8'])
% grid on
% xlabel('Number of Array Elements')
% ylabel('Simultaneous Beams')
% 
% zlabel('Power Consumption (W)')
% legend('BB Precoding','SerDes','DAC','Mixer','VCO','PS',...
%     'RF Amp','PA','Orientation','horizontal')