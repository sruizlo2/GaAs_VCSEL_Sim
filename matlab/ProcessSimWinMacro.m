function ProcessSimWinMacro(macroData, varsToPlotX, varsToPlotY,...
  plotTittle, xAxesLabel, yAxesLabel, xLims, yLims, curFig)

for thisVarX = 1:length(varsToPlotX)
  for thisVarY = 1:length(varsToPlotY{thisVarX})
    if thisVarY == 1 && length(varsToPlotY{thisVarX}) > 1
      yyaxis left
    elseif thisVarY > 1
      yyaxis right
    end
    figure(curFig + thisVarX),
    if iscell(varsToPlotY{thisVarX})
      plot(macroData.(varsToPlotX{thisVarX}), macroData.(varsToPlotY{thisVarX}{thisVarY}),...
        'linewidth', 1.5)
      ylabel(yAxesLabel{thisVarX}{thisVarY}), ylim(yLims{thisVarX}{thisVarY})
    else
      plot(macroData.(varsToPlotX{thisVarX}), macroData.(varsToPlotY{thisVarX}),...
        'linewidth', 1.5)
      ylabel(yAxesLabel{thisVarX}),  ylim(yLims{thisVarX})
    end
  end
  xlabel(xAxesLabel{thisVarX}), title(plotTittle{thisVarX})
  xlim(xLims{thisVarX}), grid on, grid minor
end
