% close all
% Pre-defined parameters for plotting
set(groot,'defaultAxesFontSize', 14)
set(groot,'defaulttextInterpreter', 'latex')
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(0,'defaultfigurecolor',[1 1 1])

%%
devIndx = 2;
if devIndx == 1
  %   filename = 'C:\SimWindows\VCSEL\StandardVCSEL_MACRO_Output.dat';
  filename = 'C:\SimWindows\VCSEL\ImprovedVCSEL_MACRO_Output.dat';
else
%   filenames = 'C:\SimWindows\VCSEL\VCSEL_MACRO_Output.dat';
  filenames = {'C:\SimWindows\VCSEL\ImprovedVCSEL_MACRO_Output.dat',...
    'C:\SimWindows\VCSEL\VCSEL_Opt_MACRO_Output.dat',...
    'C:\SimWindows\VCSEL\VCSEL_Mod1_MACRO_Output'};
end
for filename = filenames
  warning off
  macroData = readtable(filename{:}, 'VariableNamingRule', 'modify');
  warning on
  varsNames = macroData.Properties.VariableNames;
  macroData.ContJTotal_A_cm2_AtLtCont = macroData.ContJTotal_A_cm2_AtLtCont;
  macroData.ContITotal_mA_AtLtCont = macroData.ContJTotal_A_cm2_AtLtCont * 38.4e-8 * 1e6;

%   fprintf('Variaables:\n')
%   fprintf('\t%s\n', varsNames{:})

  varsToPlotX = {'ContJTotal_A_cm2_AtLtCont'};
  varsToPlotY = {{'OptPwr_mW_AtLtMir'}};
  plotTittle = {'', ''};
  xAxesLabel = {'Current density [A/cm$^2$]',...
    'Current density [A/cm$^2$]'};
  yAxesLabel ={{'Optical power [mW]'}};
  xLims = {[1e-3 1e5]};
  yLims = {[1e-12 50]};
  logScale = [true];
  curFig = 0; figure(curFig + 1), subplot(1,2,1),
  ProcessSimWinMacro(macroData, varsToPlotX, varsToPlotY,...
    plotTittle, xAxesLabel, yAxesLabel, xLims, yLims, logScale, curFig);
  hold on
  logScale = [false];
  subplot(1,2,2),
  ProcessSimWinMacro(macroData, varsToPlotX, varsToPlotY,...
    plotTittle, xAxesLabel, yAxesLabel, xLims, yLims, logScale, curFig);
  hold on
  fprintf('Optical power: %.2f mA at %.1f V\n', ...
    macroData.OptPwr_mW_AtLtMir(end), macroData.ContPot_V_AtLtCont(end))
end
subplot(1,2,1), hold off
legend('Improved VCSEL Design', 'Our Design', 'Location','northwest')
subplot(1,2,2), hold off
legend('Org.', 'conc.', 'remov. high p-doping', 'conc + remov. high p-doping', ...
  'conc + remov. high p-doping + remov. doping confinement',...
  'conc + remov. high p-doping + remov. doping confinement + lower R',...
  'conc + remov. high p-doping + remov. doping confinement + QW refractive index',...
  'conc + remov. high p-doping + remov. doping confinement + QW refractive index +  New R',...
  'Location','northwest')
