close all
% Pre-defined parameters for plotting
set(groot,'defaultAxesFontSize', 14)
set(groot,'defaulttextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])

%%
filename = 'C:\SimWindows\VCSEL\VCSEL_MACRO_Output.dat';
warning off
macroData = readtable(filename, 'VariableNamingRule', 'modify');
warning on
varsNames = macroData.Properties.VariableNames;
macroData.ContJTotal_A_cm2_AtLtCont = macroData.ContJTotal_A_cm2_AtLtCont;
macroData.ContITotal_mA_AtLtCont = macroData.ContJTotal_A_cm2_AtLtCont * 38.4e-8 * 1e6;

fprintf('Variaables:\n')
fprintf('\t%s\n', varsNames{:})

%%
varsToPlotX = {'ContJTotal_A_cm2_AtLtCont',...
               'ContJTotal_A_cm2_AtLtCont'};
varsToPlotY = {{'ContPot_V_AtLtCont', 'OptPwr_mW_AtLtMir'},...
               {'OptPwr_mW_AtLtMir'}};
plotTittle = {'', ''};
xAxesLabel = {'Current density [A/cm$^2$]',...
              'Current density [A/cm$^2$]'};
yAxesLabel ={{'Voltage [V]', 'Optical power [mW]'},...
             {'Optical power [mW]'}};
xLims = {[-inf inf],...
         [-inf inf]};
yLims = {{[-inf inf], [-inf, inf]},...
         [-inf inf]};
logScale = [false true];
curFig = 0; 
ProcessSimWinMacro(macroData, varsToPlotX, varsToPlotY,...
  plotTittle, xAxesLabel, yAxesLabel, xLims, yLims, logScale, curFig);



varsToPlotX = {'ContITotal_mA_AtLtCont',...
               'ContITotal_mA_AtLtCont'};
varsToPlotY = {{'ContPot_V_AtLtCont', 'OptPwr_mW_AtLtMir'},...
               {'OptPwr_mW_AtLtMir'}};
plotTittle = {'', ''};
xAxesLabel = {'Current [mA]',...
              'Current [mA]'};
yAxesLabel ={{'Voltage [V]', 'Optical power [mW]'},...
             {'Optical power [mW]'}};
xLims = {[-inf inf],...
         [-inf inf]};
yLims = {{[-inf inf], [-inf, inf]},...
         [-inf inf]};
logScale = [false true];
curFig = 10; 
ProcessSimWinMacro(macroData, varsToPlotX, varsToPlotY,...
  plotTittle, xAxesLabel, yAxesLabel, xLims, yLims, logScale, curFig);

