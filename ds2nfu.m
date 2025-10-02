function [xn, yn] = ds2nfu(x, y)
    ax = gca;
    axUnits = get(ax, 'Units');
    set(ax, 'Units', 'normalized');
    axPos = get(ax, 'Position');
    set(ax, 'Units', axUnits);

    xlim = get(ax, 'XLim');
    ylim = get(ax, 'YLim');

    xn = axPos(1) + (x - xlim(1)) / diff(xlim) * axPos(3);
    yn = axPos(2) + (y - ylim(1)) / diff(ylim) * axPos(4);
end
