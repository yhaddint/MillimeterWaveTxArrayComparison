% This script is used to generate fig. 8 (block power breakdown) of the paper
% Exact data is from excel sheets.

clear;clc;

DSP_FOM_range = [0.1,0.2,0.4,0.8,1.6,3.2];
font_size = 14;

% color coding in figures
grey = [0.83 0.82 0.78];
blue =  [0    0.4470    0.7410];
red = [0.8500    0.3250    0.0980];
% yellow = [ 0.9290    0.6940    0.1250];
green = [ 0    0.5    0];


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
marginw = 0.06;
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
DA_xdata_C1U8 = 64:64:256;
DA_data_C1U8 = ...
[502.2	1004.3	1506.5	2008.6;
735.8	654.7	610.0	579.4;
1264.8	2132.5	2949.1	3739.5;
640.0	1280.0	1920.0	2560.0;
3840.0	7680.0	11520.0	15360.0;
0.0	0.0	0.0	0.0;
2560.0	5120.0	7680.0	10240.0;
68996.7	32198.0	20884.7	15446.4];

DA_xdata_C1U16 = 32,64:64:256;
DA_data_C1U16 = ...
[351.5	703.0	1406.0	2109.0	2812.1;
1523.8	1335.2	1164.3	1088.0	1088.0;
531.4	901.5	1618.3	2337.8	3117.1;
320.0	640.0	1280.0	1920.0	2560.0;
1920.0	3840.0	7680.0	11520.0	15360.0;
0.0	0.0	0.0	0.0	0.0;
1280.0	2560.0	5120.0	7680.0	10240.0;
31683.1	12157.1	5103.0	3212.4	2343.3];

DA_xdata_C1U32 = 32,64:64:256;
DA_data_C1U32 = ...
[552.4	1104.7	2209.5	3314.2	4419.0;
2791.0	2353.1	2176.0	2176.0	2176.0;
472.5	814.5	1558.5	2337.8	3117.1;
320.0	640.0	1280.0	1920.0	2560.0;
1920.0	3840.0	7680.0	11520.0	15360.0;
0.0	0.0	0.0	0.0	0.0;
1280.0	2560.0	5120.0	7680.0	10240.0;
16513.1	5430.4	2258.5	1402.3	1018.2];

SA_xdata = [256 512 768 1024];
SA_data_C1U8 = ...
[175.8	326.4	477.0	627.7
1091.4	1074.8	1066.8	1062.3
1453.9	2723.0	3960.3	5188.7
320.0	640.0	960.0	1280.0
1920.0	3840.0	5760.0	7680.0
2560.0	5120.0	7680.0	10240.0
12800.0	25600.0	38400.0	51200.0
225231.1	95066.7	58455.9	41862.8];

SA_data_C1U16 = ...
[251.1	401.7	552.4	703.0
1809.0	1702.4	1667.6	1652.3
757.2	1306.6	1874.9	2453.2
320.0	640.0	960.0	1280.0
1920.0	3840.0	5760.0	7680.0
2560.0	5120.0	7680.0	10240.0
12800.0	25600.0	38400.0	51200.0
134843.0	39249.0	21919.4	15210.0];

SA_xdata_U32 = [64,128,256,512];
SA_data_C1U32 = ...
[439.4	477.0	552.4	703.0
2176.0	2176.0	2351.3	2552.6
97.4	194.8	407.1	865.0
80.0	160.0	320.0	640.0
480.0	960.0	1920.0	3840.0
640.0	1280.0	2560.0	5120.0
3200.0	6400.0	12800.0	25600.0
28302.7	26596.7	21618.6	18023.0];


FH_xdata = 64:64:256;
FH_data_C1U8 = ...
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


SA_xdata_C2U2 = [256 512 768 1024];
SA_data_C2U2 = ...
[453.5	604.2	905.4	1206.7;
171.8	172.2	172.1	172.1;
1393.2	1862.2	2791.4	3721.8;
960.0	1280.0	1920.0	2560.0;
5760.0	7680.0	11520.0	15360.0;
7680.0	10240.0	15360.0	20480.0;
38400.0	51200.0	76800.0	102400.0;
168588.6	128477.9	85276.3	63948.2];

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
[0	0	0	0
0.0	0.0	0.0	0.0
143.0	143.0	143.2	143.2
10.0	10.0	10.0	10.0
60.0	60.0	60.0	60.0
3840.0	5120.0	7680.0	10240.0
25600.0	34133.3	51200.0	68266.7
68364.1	51265.9	34263.2	25693.8];
%%
for ii=1:length(DSP_FOM_range)
    if ii>1
        DSP_scaling_old = DSP_FOM_range(ii-1);
        DSP_scaling_new = DSP_FOM_range(ii);;
    else
        DSP_scaling_old = 1/13;
        DSP_scaling_new = DSP_FOM_range(ii);
    end

    DA_data_C1U8(1,:) = DA_data_C1U8(1,:)*(DSP_scaling_new/DSP_scaling_old);
    DA_data_C1U16(1,:) = DA_data_C1U16(1,:)*(DSP_scaling_new/DSP_scaling_old);
    DA_data_C1U32(1,:) = DA_data_C1U32(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    DA_data_C2U2(1,:) = DA_data_C2U2(1,:)*(DSP_scaling_new/DSP_scaling_old);
    DA_data_C2U4(1,:) = DA_data_C2U4(1,:)*(DSP_scaling_new/DSP_scaling_old);
    DA_data_C2U8(1,:) = DA_data_C2U8(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    DA_data_C3U1(1,:) = DA_data_C3U1(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SA_data_C1U8(1,:) = SA_data_C1U8(1,:)*(DSP_scaling_new/DSP_scaling_old);
    SA_data_C1U16(1,:) = SA_data_C1U16(1,:)*(DSP_scaling_new/DSP_scaling_old);
    SA_data_C1U32(1,:) = SA_data_C1U32(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    SA_data_C2U2(1,:) = SA_data_C2U2(1,:)*(DSP_scaling_new/DSP_scaling_old);
    SA_data_C2U4(1,:) = SA_data_C2U4(1,:)*(DSP_scaling_new/DSP_scaling_old);
    SA_data_C2U8(1,:) = SA_data_C2U8(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    SA_data_C3U1(1,:) = SA_data_C3U1(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FH_data_C1U8(1,:) = FH_data_C1U8(1,:)*(DSP_scaling_new/DSP_scaling_old);
    FH_data_C1U16(1,:) = FH_data_C1U16(1,:)*(DSP_scaling_new/DSP_scaling_old);
    FH_data_C1U32(1,:) = FH_data_C1U32(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    FH_data_C2U2(1,:) = FH_data_C2U2(1,:)*(DSP_scaling_new/DSP_scaling_old);
    FH_data_C2U4(1,:) = FH_data_C2U4(1,:)*(DSP_scaling_new/DSP_scaling_old);
    FH_data_C2U8(1,:) = FH_data_C2U8(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    FH_data_C3U1(1,:) = FH_data_C3U1(1,:)*(DSP_scaling_new/DSP_scaling_old);
    
    [DA_bestC1(ii), DA_indexC1(ii)] = min([sum(DA_data_C1U8,1),sum(DA_data_C1U16,1),sum(DA_data_C1U32,1)]);
    [SA_bestC1(ii), SA_indexC1(ii)] = min([sum(SA_data_C1U8,1),sum(SA_data_C1U16,1),sum(SA_data_C1U32(:,3:4),1)]);
    [FH_bestC1(ii), FH_indexC1(ii)] = min([sum(FH_data_C1U8,1),sum(FH_data_C1U16,1),sum(FH_data_C1U32,1)]);
    
    [DA_bestC2(ii), DA_indexC2(ii)] = min([sum(DA_data_C2U2,1),sum(DA_data_C2U4,1),sum(DA_data_C2U8,1)]);
    [SA_bestC2(ii), SA_indexC2(ii)] = min([sum(SA_data_C2U2,1),sum(SA_data_C2U4,1),sum(SA_data_C2U8,1)]);
    [FH_bestC2(ii), FH_indexC2(ii)] = min([sum(FH_data_C2U2,1),sum(FH_data_C2U4,1),sum(FH_data_C2U8,1)]);
    
    [DA_bestC3(ii), DA_indexC3(ii)] = min([sum(DA_data_C3U1,1)]);
    [SA_bestC3(ii), SA_indexC3(ii)] = min([sum(SA_data_C3U1,1)]);
    [FH_bestC3(ii), FH_indexC3(ii)] = min([sum(FH_data_C3U1,1)]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   case I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
h= subplot('position', [axisl(1), axisb(1), axisw, axish] );
gcf = semilogx(DSP_FOM_range*10,DA_bestC1/1e3,'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = semilogx(DSP_FOM_range*10,SA_bestC1/1e3,'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = semilogx(DSP_FOM_range*10,FH_bestC1/1e3,'-x','linewidth',2,'markersize',8);hold on
gcf.Color = green;
grid on
xticks([1,2,4,8,16,32])
yticks([0,25,50,75])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',['0.1';'0.2';'0.4';'0.8';'1.6';'3.2'])
xlim([1/sqrt(2),32*sqrt(2)])   
ylim([0,80])                  
xlabel('DSP Energy Efficiency (mW/GOPS)')
ylabel('Power Consumption (W)')

for ii=1:length(DSP_FOM_range)
    
    % ------------- DA ---------------
    if ii<5
        txtDA = '$\{64,32\}$';
    else
        txtDA = '$\{32,32\}$';
    end
    
    text(DSP_FOM_range(ii)*10,DA_bestC1(ii)/1e3-7.5,txtDA,...
        'HorizontalAlignment','center','FontSize', font_size,'Color',red)
    
    % ------------- SA ---------------
    txtSA = '$\{256,32\}$';
    if ii<5
    text(DSP_FOM_range(ii)*10,SA_bestC1(ii)/1e3+5,txtSA,...
        'HorizontalAlignment','center','FontSize', font_size,'Color',blue)
    else
        text(DSP_FOM_range(ii)*10,SA_bestC1(ii)/1e3-7.5,txtSA,...
        'HorizontalAlignment','center','FontSize', font_size,'Color',blue)
    end
    
    % ------------- FH ---------------
    txtFH = '$\{64,16\}$';
    text(DSP_FOM_range(ii)*10,FH_bestC1(ii)/1e3+5,txtFH,...
        'HorizontalAlignment','center','FontSize', font_size,'Color',green)
end
t = title('Dense Urban MBB',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   case II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h= subplot('position', [axisl(2), axisb(2), axisw, axish] );
gcf = semilogx(DSP_FOM_range*10,DA_bestC2/1e3,'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = semilogx(DSP_FOM_range*10,SA_bestC2/1e3,'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = semilogx(DSP_FOM_range*10,FH_bestC2/1e3,'-x','linewidth',2,'markersize',8);hold on
gcf.Color = green;
grid on
xticks([1,2,4,8,16,32])
yticks([100,125,150,175,200,225,250])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',['0.1';'0.2';'0.4';'0.8';'1.6';'3.2'])
ylim([100,250])
xlim([1/sqrt(2),55]) 
xlabel('DSP Energy Efficiency (mW/GOPS)')
ylabel('Power Consumption (W)')
for ii=1:length(DSP_FOM_range)
    
    % ------------- DA ---------------
    if ii<6
        txtDA = '$\{384,8\}$';
    else
        txtDA = '$\{256,8\}$';
    end
    if ii==5
        text(20.35,163,txtDA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',red)
    elseif ii==6
        text(35.2,218.2,txtDA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',red)
    else
        text(DSP_FOM_range(ii)*10,DA_bestC2(ii)/1e3-10,txtDA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',red)
    end
    
    % ------------- SA ---------------
    txtSA = '$\{1536,2\}$';
    if ii<7
        text(DSP_FOM_range(ii)*10,SA_bestC2(ii)/1e3+10,txtSA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',blue) 
    else
        text(37.1,230,txtSA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',blue)
    end
    
    % ------------- FH ---------------
    txtFH = '$\{768,2\}$';
    if ii==5    
        text(13,180,txtFH,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',green)  
    elseif ii==6
        text(32,180,txtFH,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',green)
    else
        text(DSP_FOM_range(ii)*10,FH_bestC2(ii)/1e3-10,txtFH,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',green)
    end

end
t = title('50+Mbps Everywhere',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   case III
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h= subplot('position', [axisl(3), axisb(3), axisw, axish] );
gcf = semilogx(DSP_FOM_range*10,DA_bestC3/1e3,'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = semilogx(DSP_FOM_range*10,SA_bestC3/1e3,'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = semilogx(DSP_FOM_range*10,FH_bestC3/1e3,'-x','linewidth',2,'markersize',8);hold on
gcf.Color = green;
grid on
xticks([1,2,4,8,16,32])
yticks([75,100,125,150,175])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',['0.1';'0.2';'0.4';'0.8';'1.6';'3.2'])

ylim([75,175])
xlim([1/sqrt(2),32*sqrt(2)]) 
xlabel('DSP Energy Efficiency (mW/GOPS)')
ylabel('Power Consumption (W)')
for ii=1:length(DSP_FOM_range)

    % ------------- DA ---------------
    if ii<6
        txtDA = '$\{384,1\}$';
        text(DSP_FOM_range(ii)*10,DA_bestC3(ii)/1e3+7.5,txtDA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',red)
    else
        txtDA = '$\{256,1\}$';
        text(20,170,txtDA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',red)
    end

    % ------------- SA ---------------
    txtSA = '$\{512,1\}$';
    if ii<6
        text(DSP_FOM_range(ii)*10,SA_bestC3(ii)/1e3+5,txtSA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',blue)
    else
        text(DSP_FOM_range(ii)*10,SA_bestC3(ii)/1e3+5,txtSA,...
            'HorizontalAlignment','center','FontSize', font_size,'Color',blue)
    end
    
    % ------ FH ------------
    text(DSP_FOM_range(ii)*10,FH_bestC3(ii)/1e3-5,txtSA,...
        'HorizontalAlignment','center','FontSize', font_size,'Color',green)
    
end
t = title('Self-backhauling',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position
% ------------- Legends --------------
legend('DA','SA','FH','Orientation','horizontal')
set(legend,'Position',[0.4, axish+0.175, axisw, 0.05])


