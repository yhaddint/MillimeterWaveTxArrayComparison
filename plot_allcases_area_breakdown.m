clear;clc;
%%
grey = [0.83 0.82 0.78];
green = [0.47 0.67 0.19];
scarlet = [0.64 0.08 0.18];
carbon = [0.5,0.5,0.5];
pool = [0.3,0.75,0.93];

%%
DA_xdata = [128,256,384];
DA_data = ...
[36.6	73.1	109.7;
2.5	2.5	2.5;
6.4	12.8	19.2;
23.0	46.1	69.1;
0.0	0.0	0.0;
4.8	9.9	15.0;
3.2	6.4	9.6];

SA_xdata = 256:256:1024;
SA_data = ...
[14.4	28.7	43.1	57.4;
2.5	2.5	2.5	2.5;
1.6	3.2	4.8	6.4;
5.8	11.5	17.3	23.0;
20.5	41.0	61.4	81.9;
10.2	20.4	30.7	40.9;
8.0	16.0	24.0	32.0];

FH_xdata = [128,256,384];
FH_data = ...
[4.6	4.6	4.6	4.6;
0.0	0.0	0.0	0.0;
0.4	0.4	0.4	0.4;
1.4	1.4	1.4	1.4;
81.9	163.8	245.8	327.7;
76.5	153.3	230.1	306.9;
20.3	40.5	60.8	81.1];
%%
figure(1)

text_loc = 700;
font_size = 14;
figure
% ------------- Subplot for all Arrays --------------
b = bar(1:13,[DA_data(:,1:3).';zeros(1,7);SA_data(:,1:4).';zeros(1,7);FH_data(:,1:4).'],'stacked');
colormap parula
b(4).FaceColor = green;
b(5).FaceColor = pool;
b(6).FaceColor = carbon;
b(7).FaceColor = scarlet;
xticks([1,2,3,5,6,7,8,10,11,12,13])
set(gca,'FontSize', font_size)
set(gca,'xticklabel',[' 128';' 256';' 384';' 256';' 512';' 768';'1024';' 128';' 256';' 384';' 512'])
xlim([0.5,13.5])
ylim([0,1000])
grid on
xlabel('Number of Array Elements')
ylabel('Estimated Total IC Area (mm$^2$)')
% t = title('Digital Array',...
%     'Units', 'normalized',...
%     'Position', [0.5, 0.9, 0]); % Set Title with correct Position

text(2,text_loc,'DA','HorizontalAlignment','center','Fontsize',font_size)
text(6.5,text_loc,'SA','HorizontalAlignment','center','Fontsize',font_size)
text(11.5,text_loc,'FH','HorizontalAlignment','center','Fontsize',font_size)

% ------------- Legends --------------
legend('BB Precoding','SerDes','DAC','Mixer & VCO','PS','Split. & Comb.',...
    'RF Amp','Orientation','horizontal')
set(legend,'Position',[0.4, 0.85, 0.25, 0.05],'NumColumns',4)
