% This script is used to plot figures for required transmission power
% The data is from simulation (xx_power_theo.m)


clear;clc;clf

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

%% -------------------------------------------------------------------
%                          Case I
% -------------------------------------------------------------------


SA = zeros(7,2);
FH = zeros(8,2);

DA_size_C1U8 = [64, 128, 192, 256, 384, 512, 768, 1024]';
DA_pow_C1U8 = 46 + [-4.9405   -8.2475  -10.1271  -11.4454  -13.2510  -14.5284  -16.3367  -17.6083]';
% From DA_required_power_theo.m with SNR_origin = 18.7, target SINR = 22.1
% such target is to meet 58.5bps/Hz SE with U = 8

DA_size_C1U16 = [32, 64, 128, 192, 256, 384, 512, 768, 1024]';
DA_alpha_C1U16 = [-5];
DA_pow_C1U16 = 46 + [-8.32, -12.4857  -16.2512  -18.2642  -19.6343  -21.5234  -22.8441  -24.6672  -25.9392]';
% From DA_required_power_theo.m with SNR_origin = 18.7, target SINR = 10.71
% such target is to meet 58.5bps/Hz SE with U = 16

DA_size_C1U32 = [32, 64, 128, 192, 256, 384, 512, 768, 1024]';
DA_alpha_C1U32 = [-5];
DA_pow_C1U32 = 46 + [-11.15, -15.9819  -19.7928  -21.8601  -23.2542  -25.1158  -26.4309  -28.2500  -29.5333]';
% From DA_required_power_theo.m with SNR_origin = 18.7, target SINR = 4.11
% such target is to meet 58.5bps/Hz SE with U = 32

SA_C1U8(:,1) = [64, 128, 256, 512, 768, 1024]';
SA_C1U8(:,2) = 46 + [8.6833    4.5097    0.1984   -3.5480   -5.6640   -7.1179]';
% From SA_required_power_theo.m with SNR_origin = 18.7, target SINR = 22.1
% such target is to meet 58.5bps/Hz SE with U = 8

SA_C1U16(:,1) = [ 64, 128, 256, 512, 768, 1024]';
SA_C1U16(:,2) = 46 + [2.3169    0.3376   -2.7350   -7.3995   -9.9247  -11.5067]';
% From SA_required_power_theo.m with SNR_origin = 18.7, target SINR = 10.1
% such target is to meet 58.5bps/Hz SE with U = 16

SA_C1U32(:,1) = [ 64, 128, 256, 512, 768, 1024]';
SA_C1U32(:,2) = 46 + [-8.8129   -9.0949   -9.8275  -10.7768  -11.8748  -12.8552]';
% From SA_required_power_theo.m with SNR_origin = 18.7, target SINR = 4
% such target is to meet 58.5bps/Hz SE with U = 32

FH_size_C1U8 = [64, 128, 192, 256, 384, 512, 768, 1024]';
FH_alpha_C1U8 = []';
FH_pow_C1U8 = 46 + [-4.1682   -7.4068   -9.2433  -10.5359  -12.3323  -13.5959  -15.3800  -16.6388 ]';
% From FH_required_power_theo.m with SNR_origin = 18.7, target SINR = 22.1
% such target is to meet 58.5bps/Hz SE with U = 8

FH_size_C1U16 = [32, 64, 128, 192, 256, 384, 512, 768, 1024]';
FH_alpha_C1U16 = [0     8    16    20    22    26    28    30    30]';
FH_pow_C1U16 = 46 + [-6.6289, -11.3315  -15.2515  -17.2832  -18.6557  -20.5481  -21.8811  -23.6826  -24.9589]';
% From FH_required_power_theo.m with SNR_origin = 18.7, target SINR = 22.1
% such target is to meet 58.5bps/Hz SE with U = 16

FH_size_C1U32 = [32, 64, 128, 192, 256, 384, 512, 768, 1024]';
FH_alpha_C1U32 = [2    12    20    24    26    30    30    30    30]';
FH_pow_C1U32 = 46 + [-8.5469  -13.9579 -18.2311  -20.4983  -21.9683  -23.9177  -25.2801  -27.1228  -28.4239]';
% From FH_required_power_theo.m with SNR_origin = 18.7, target SINR = 22.1
% such target is to meet 58.5bps/Hz SE with U = 16

figure(1)
h= subplot('position', [axisl(1), axisb(1), axisw, axish] );
set(0,'defaultTextInterpreter','latex'); %trying to set the default
% semilogx(DA(:,1),DA(:,2),'-o','linewidth',2);hold on
% semilogx(SA(:,1),SA(:,2),'-o','linewidth',2);hold on
% semilogx(FH(:,1),FH(:,2),'-o','linewidth',2);hold on
gcf = plot(DA_size_C1U8,DA_pow_C1U8,'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(DA_size_C1U16,DA_pow_C1U16,'--o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(DA_size_C1U32,DA_pow_C1U32,'-.o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(SA_C1U8(:,1),SA_C1U8(:,2),'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(SA_C1U16(:,1),SA_C1U16(:,2),'--s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(SA_C1U32(:,1),SA_C1U32(:,2),'-.s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(FH_size_C1U8,FH_pow_C1U8,'-x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(FH_size_C1U16,FH_pow_C1U16,'--x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(FH_size_C1U32,FH_pow_C1U32,'-.x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
grid on

set(gca,'xtick',[64,128:128:1024])
set(gca,'FontSize', 14)

xlim([64,1024])
% ylim([30,55])
     
xlabel('Number of Array Elements')
ylabel('Required Transmission Power (dBm)')
% title('Dense Urban (LOS 100m, 58.8bps/Hz)')
t = title('Dense Urban (LOS 100m, 58.8bps/Hz)',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

legend('DA (U = 8)','DA (U = 16)','DA (U = 32)',...
       'SA (U = 8)', 'SA (U = 16)','SA (U = 32)',...
       'FH (U = 8)','FH (U = 16)','FH (U = 32)')
set(legend,'NumColumns',3);
set(legend,'Position',[axisl(1), axish+0.15, axisw, 0.15])
%% -------------------------------------------------------------------
%                          Case II
% -------------------------------------------------------------------
DA = zeros(8,4);
SA = zeros(8,2);
FH = zeros(8,2);
DA(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';

DA(:,2) = 46 + [6.0134    2.9005    1.0771   -0.1957   -1.9656   -3.2242   -4.9993, -6.19]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = 6.2
% such target is to meet 4.7bps/Hz SE with U = 2

DA(:,3) = 46 + [3.8027    0.7891   -0.9719   -2.2548   -4.0183   -5.2811   -7.0536, -8.39]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = 1.0
% such target is to meet 4.7bps/Hz SE with U = 4

DA(:,4) = 46 + [2.3864   -0.5175   -2.2236   -3.4484   -5.1747   -6.3939   -8.1201, -9.39]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = -3.0
% such target is to meet 4.7bps/Hz SE with U = 8

SA_C2U2(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048]';
SA_C2U2(:,2) = 46 + [9.7333    6.7747    4.9415    3.5939    1.8716    0.5997   -1.0659   -2.2412, -4.0252   -5.2730]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = 6.2
% such target is to meet 4.7bps/Hz SE with U = 2

SA_C2U4(:,1) = [64, 128, 256, 384, 512, 768, 1024, 1536, 2048]';
SA_C2U4(:,2) = 46 + [10.7234    7.7125    4.6322    2.8522    1.5527   -0.1843   -1.3458  -3.0708   -4.3298]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = 1.0
% such target is to meet 4.7bps/Hz SE with U = 4

SA_C2U8(:,1) = [64, 128, 256, 512, 768, 1024, 1536, 2048]';
SA_C2U8(:,2) = 46 + [10.0185    8.6198    6.0829    3.4321    1.4375    0.3329  -1.2521   -2.3857]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = -3.0
% such target is to meet 4.7bps/Hz SE with U = 8

FH(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';
FH(:,2) = 46 + [6.6914    3.7013    1.9801    0.7487   -1.0217   -2.2829   -4.0169   -5.2696]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = 6.2
% such target is to meet 4.7bps/Hz SE with U = 2

FH(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';
FH(:,3) = 46 + [4.8377    1.7614   -0.0189   -1.3021   -3.0557   -4.3175   -6.0891   -7.3483]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = 1.0
% such target is to meet 4.7bps/Hz SE with U = 4

FH(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';
FH(:,4) = 46 + [3.6985    0.7333   -1.0634   -2.3021   -4.0412   -5.2841   -7.0515   -8.3137]';
% From DA_required_power_theo.m with SNR_origin = -14.7, target SINR = -3.0
% such target is to meet 58.5bps/Hz SE with U = 8


h= subplot('position', [axisl(2), axisb(2), axisw, axish] );
set(0,'defaultTextInterpreter','latex'); %trying to set the default
% semilogx(DA(:,1),DA(:,2),'-o','linewidth',2);hold on
% semilogx(SA(:,1),SA(:,2),'-o','linewidth',2);hold on
% semilogx(FH(:,1),FH(:,2),'-o','linewidth',2);hold on
gcf = plot(DA(:,1),DA(:,2),'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(DA(:,1),DA(:,3),'--o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(DA(:,1),DA(:,4),'-.o','linewidth',2,'markersize',8);hold on
gcf.Color = red;

gcf = plot(SA_C2U2(:,1),SA_C2U2(:,2),'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(SA_C2U4(:,1),SA_C2U4(:,2),'--s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
gcf = plot(SA_C2U8(:,1),SA_C2U8(:,2),'-.s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;

gcf = plot(FH(:,1),FH(:,2),'-x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(FH(:,1),FH(:,3),'--x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;
gcf = plot(FH(:,1),FH(:,4),'-.x','linewidth',2,'markersize',8);hold on
gcf.Color = yellow;

grid on
set(gca,'xtick',[64,128:128:1024])
set(gca,'FontSize', 14)

xlabel('Number of Array Elements')
ylabel('Required Transmission Power (dBm)')
xlim([64,1024])
t = title('50+Mbps Everywhere (NLOS 100m, 4.7bps/Hz)',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position

% ylim([30,55])
% legend('DA', 'SA', 'FH')
legend('DA (U = 2)','DA (U = 4)','DA (U = 8)',...
       'SA (U = 2)', 'SA (U = 4)', 'SA (U = 8)',...
       'FH (U = 2)','FH (U = 4)','FH (U = 8)')
set(legend,'NumColumns',3);
set(legend,'Position',[axisl(2), axish+0.15, axisw, 0.15])
%% -------------------------------------------------------------------
%                          Case III
% -------------------------------------------------------------------
DA = zeros(8,2);
SA = zeros(8,2);
FH = zeros(8,2);
DA(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';
SA(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';
FH(:,1) = [64, 128, 192, 256, 384, 512, 768, 1024]';

DA(:,2) = 46 + [1.8422   -1.1681   -2.9290   -4.1784   -5.9393   -7.1887   -8.9496  -10.1990]';
SA(:,2) = 46 + [2.5958   -0.2709   -1.9815   -3.2461   -4.9825   -6.2258   -7.9829   -9.2266]';


h= subplot('position', [axisl(3), axisb(3), axisw, axish] );
set(0,'defaultTextInterpreter','latex'); %trying to set the default
% semilogx(DA(:,1),DA(:,2),'-o','linewidth',2);hold on
% semilogx(SA(:,1),SA(:,2),'-o','linewidth',2);hold on
% semilogx(FH(:,1),FH(:,2),'-o','linewidth',2);hold on
gcf = plot(DA(:,1),DA(:,2),'-o','linewidth',2,'markersize',8);hold on
gcf.Color = red;
gcf = plot(SA(:,1),SA(:,2),'-s','linewidth',2,'markersize',8);hold on
gcf.Color = blue;
grid on
set(gca,'xtick',[64,128:128:1024])
set(gca,'FontSize', 14)
xlabel('Number of Array Elements')
ylabel('Required Transmission Power (dBm)')
xlim([64,1024])
t = title('Self-Backhauling (LOS 707m, 11.8bps/Hz)',...
    'Units', 'normalized',...
    'Position', [0.5, 0.9, 0]); % Set Title with correct Position


% ylim([30,55])
legend('DA (U = 1)', 'SA and FH (U = 1)')
set(legend,'NumColumns',3);
set(legend,'Position',[axisl(3), axish+0.2, axisw, 0.05])
%%
