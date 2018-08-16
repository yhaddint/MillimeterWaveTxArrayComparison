clear;clc;
% CMOS blue
% BiCMOS red
% GaAs 
loc_os_x = 0.5;
loc_os_y = 0.3;
label = cell(12,9);


label(1,:) = {15.1, 33.7, '[69]', '2017','CMOS','40nm','Shakib', -5, -2};
label(2,:) = {18.7, 12.4, '[70]', '2017','CMOS','28nm','Moret', loc_os_x, loc_os_y};
label(3,:) = {14, 35.5, '[71]', '2016','CMOS','28nm','Shakib', -3.5, 2};
label(4,:) = {19.8, 43.3, '[72]', '2016','CMOS','28nm','Park', loc_os_x, loc_os_y};

label(5,:) = {18.6, 35.3, '[73]', '2017','BiCMOS','0.13$\mu$m','Sarkar', loc_os_x, loc_os_y};
label(6,:) = {23.7, 32, '[74]', '2017','BiCMOS','0.18$\mu$m','Sarkar', loc_os_x, loc_os_y};
label(7,:) = {20.6, 29, '[75]', '2014','BiCMOS','0.18$\mu$m','', -3.5, -2};

label(8,:) = {26, 40, '[76]', '2017','GaAs','','Nguyen', loc_os_x, loc_os_y};
label(9,:) = {32, 30, '[77]', 'TGA4544','GaAs','','Qorvo', -1, -1.5};
label(10,:) = {30.5, 22, '[78]', 'HMC1132','GaAs','','AD', -4, -2};

label(11,:) = {38.5, 24, '[79]', '2015','GaN','','', -1, -2};
label(12,:) = {45.5, 32, '[80]', '2015','GaN','','', -3, -2};


figure
for ii=1:12
    switch label{ii,5}
        case 'CMOS'
            setcolor = [0    0.4470    0.7410];
        case 'BiCMOS'
            setcolor = [0.8500    0.3250    0.0980];
        case 'GaAs'
            setcolor = [0.4660    0.6740    0.1880];
        case 'GaN'
            setcolor = 'k';
    end
p(ii) = plot(label{ii,1},label{ii,2},'o','MarkerSize',8,'MarkerEdgeColor',setcolor, 'MarkerFaceColor',setcolor);
hold on
t = text(label{ii,1}+label{ii,8},label{ii,2}+label{ii,9},[label{ii,3} ', ' label{ii,4} ', ' label{ii,6}],'Interpreter','latex');
t.FontSize = 14;
t.Color = setcolor;
end

xlabel('Saturated Output Power (dBm)')
ylabel('Power-Adding Efficiency (%)')
grid on
legend([p(1),p(5),p(8),p(11)],'CMOS','SiGe BiCMOS','GaAs','GaN');