function ProcessSimWinMacro(macroData, varsToPlotX, varsToPlotY,...
  plotTittle, xAxesLabel, yAxesLabel, xLims, yLims, logScale, curFig)

for thisVarX = 1:length(varsToPlotX)
  figure(curFig + thisVarX),
  for thisVarY = 1:length(varsToPlotY{thisVarX})
    if thisVarY == 1 && length(varsToPlotY{thisVarX}) > 1
      yyaxis left
    elseif thisVarY > 1
      yyaxis right
    end
    if logScale(thisVarX)
      loglog(macroData.(varsToPlotX{thisVarX}), macroData.(varsToPlotY{thisVarX}{thisVarY}),...
        'linewidth', 1.5)
    else
      plot(macroData.(varsToPlotX{thisVarX}), macroData.(varsToPlotY{thisVarX}{thisVarY}),...
        'linewidth', 1.5)
    end
    if length(varsToPlotY{thisVarX}) > 1
      ylabel(yAxesLabel{thisVarX}{thisVarY}), ylim(yLims{thisVarX}{thisVarY})
    else
      ylabel(yAxesLabel{thisVarX}),  ylim(yLims{thisVarX})
    end
  end
  xlabel(xAxesLabel{thisVarX}), title(plotTittle{thisVarX})
  xlim(xLims{thisVarX}), grid on, grid minor
end
