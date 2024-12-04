% close all
% Pre-defined parameters for plotting
set(groot,'defaultAxesFontSize', 14)
set(groot,'defaulttextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])

%% 
filename = 'C:\SimWindows\VCSEL\VCSEL_OrgDesign_ModeField.dat';
warning off
macroData = readtable(filename, 'VariableNamingRule', 'modify');
warning on

figure(100), plot(macroData.ModeTotFieldMag_micr_3_2_, macroData.GridPos_microns_, 'k')
set(gca, 'YDir','reverse')
ylabel('$x$ $\mu$m'), xlabel('Field Amplitude'), grid on, grid minor
axis tight