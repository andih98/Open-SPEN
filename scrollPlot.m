function scrollPlot(data, dim, varargin)
% interactiveScrollND Interaktives Durchblättern von Schnitten in einem N-D Array.
%
%   interactiveScrollND(data, dim) öffnet eine GUI-Figur mit Slider und
%   erlaubt, durch die Indizes entlang Dimension dim zu scrollen. Bei
%   Vektoren wird ein 1D-Plot (plot) erzeugt, bei 2D-Slices imagesc. In
%   höheren Dimensionen wird squeeze verwendet, zusätzliche Dimensionen
%   müssen per Skript zuvor zusammengefasst werden.
%
%   interactiveScrollND(data, dim, 'Display','abs',...)  legt fest, wie
%   komplexe Daten angezeigt werden. Mögliche Optionen:
%       'abs'  -> Betrag  (Default bei 2D-Images)
%       'real' -> Realteil (Default bei 1D, falls nicht spezifiziert)
%       'imag' -> Imaginärteil
%       'angle'-> Phase (Winkel)
%   Weitere Name-Value-Paare:
%       'TitlePrefix' : Text, der vor "Slice k/N" im Titel angezeigt wird.
%       'Colormap'    : Z.B. 'gray', 'jet' etc. (für 2D imagesc).
%       'CLim'        : Zweiwertiger Vektor [cmin cmax] zur Farbskalierung.
%                        Wenn leer (Default), wird auto gewählt.
%
%   Beispiel:
%       vol = randn(64,64,30) + 1i*randn(64,64,30);
%       interactiveScrollND(vol, 3, 'Display','abs', 'Colormap','gray');
%
%   Eingabe:
%     data : N-dimensionales Array (beliebiger Typ: real oder komplex).
%     dim  : die Dimension, entlang derer gescrollt werden soll (1 ≤ dim ≤ ndims(data)).
%
%   Die GUI enthält:
%     - Slider für Index von 1 bis size(data,dim).
%     - Titel: zeigt "TitlePrefix Slice k/N".
%     - KeyPressFcn: linke/rechte Pfeile ändern Index um ±1.
%     - Automatische Erkennung: 1D-Slice => plot, 2D-Slice => imagesc.
%
%   Hinweis:
%   - Wenn data z.B. 4D ist, und Sie entlang dim scrollen, muss der
%     verbleibende Slice (nach squeeze) 1D oder 2D sein. Andernfalls
%     (z.B. 3D-Rest) müsste man weitere Parameter übergeben, um zu
%     entscheiden, welche Ebene man in den 3D-Slices anzeigen will.
%
%   Autor: ChatGPT-Beispiel, 2025

    % Input-Parser für Name-Value-Paare
    p = inputParser;
    validDisplay = {'abs','real','imag','angle'};
    addRequired(p, 'data');
    addRequired(p, 'dim', @(x) isnumeric(x) && isscalar(x) && x>=1 && mod(x,1)==0);
    addParameter(p, 'Display', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'TitlePrefix', '', @(s) ischar(s) || isstring(s));
    addParameter(p, 'Colormap', 'parula', @(s) ischar(s) || isstring(s) || isempty(s));
    addParameter(p, 'CLim', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
    parse(p, data, dim, varargin{:});
    data = p.Results.data;
    dim   = p.Results.dim;
    dispOpt = char(p.Results.Display);
    titlePrefix = char(p.Results.TitlePrefix);
    cmap = char(p.Results.Colormap);
    clim = p.Results.CLim;

    % Prüfen von dim
    nd = ndims(data);
    if dim > nd
        error('Dimension dim (%d) größer als ndims(data) (%d).', dim, nd);
    end
    N = size(data, dim);
    if N < 1
        error('Die gewählte Dimension hat Länge < 1.');
    end

    % GUI-Figur und Achse anlegen
    hFig = figure('Name','Interactive Scroll', ...
                  'NumberTitle','off', ...
                  'KeyPressFcn', @onKeyPress);
    hAx = axes('Parent', hFig);
    
    % Initial-Index
    currIdx = 1;
    
    % Initial-Slice extrahieren und plotten
    slice0 = getSlice(data, dim, currIdx);  % squeeze
    [plotHandle, sliceType] = plotSlice(hAx, slice0, dispOpt, cmap, clim);
    updateTitle();

    % Slider anlegen
    hSlider = uicontrol('Parent', hFig, 'Style','slider', ...
                        'Units','normalized', ...
                        'Position',[0.15 0.02 0.7 0.04], ...
                        'Min', 1, 'Max', N, 'Value', currIdx, ...
                        'SliderStep', [1/(N-1), min(10, N-1)/(N-1)], ...
                        'Callback', @onSliderChange);
    % Textfeld zur Anzeige "k/N"
    hText = uicontrol('Parent', hFig, 'Style','text', ...
                      'Units','normalized', ...
                      'Position',[0.87 0.02 0.12 0.04], ...
                      'String', sprintf('%d/%d', currIdx, N));
    
    % Callback Slider
    function onSliderChange(src, ~)
        k = round(get(src,'Value'));
        k = max(1, min(N, k));
        set(src, 'Value', k);
        updateSlice(k);
    end

    % Callback Tastatur
    function onKeyPress(~, event)
        switch event.Key
            case 'rightarrow'
                newIdx = currIdx + 1;
            case 'leftarrow'
                newIdx = currIdx - 1;
            otherwise
                return;
        end
        newIdx = max(1, min(N, newIdx));
        if newIdx ~= currIdx
            set(hSlider, 'Value', newIdx);
            updateSlice(newIdx);
        end
    end

    % Funktion: Slice aktualisieren
    function updateSlice(k)
        currIdx = k;
        slice = getSlice(data, dim, currIdx);
        % Aktualisiere Plot je nach Typ:
        switch sliceType
            case 'vector'
                ydata = mapComplex(slice, dispOpt);
                set(plotHandle, 'YData', ydata);
                % xData bleibt automatisch 1:length(slice)
                ylim(hAx, 'auto');
            case '2D'
                img = mapComplex(slice, dispOpt);
                set(plotHandle, 'CData', img);
                % evtl. Achsenlimits fixieren oder auto:
                if isempty(clim)
                    set(hAx, 'CLimMode', 'auto');
                else
                    caxis(hAx, clim);
                end
            otherwise
                % Sollte nicht vorkommen
                error('Unbekannter sliceType.');
        end
        % Titel und Text aktualisieren
        updateTitle();
        set(hText, 'String', sprintf('%d/%d', currIdx, N));
    end

    function updateTitle()
        if isempty(titlePrefix)
            title(hAx, sprintf('Slice %d / %d (dim %d)', currIdx, N, dim));
        else
            title(hAx, sprintf('%s Slice %d / %d (dim %d)', titlePrefix, currIdx, N, dim));
        end
    end
end

% Hilfsfunktion: Extrahiere die k-te Ebene entlang dimension dim und squeeze
function sl = getSlice(data, dim, k)
    % Erzeuge Index-Zell-Array: ':' für alle Dimensionen außer dim, wo k
    idx = repmat({':'}, 1, ndims(data));
    idx{dim} = k;
    sl = squeeze(data(idx{:}));
    % sl kann nun 0D (Skalar), 1D (Vektor) oder 2D (Matrix) sein
    % Skalar behandeln wir wie Vektor der Länge 1.
end

% Hilfsfunktion: Plotte slice initial, gib Handle und sliceType zurück
function [hPlot, sliceType] = plotSlice(ax, slice, dispOpt, cmap, clim)
    % Entscheide, ob vector oder 2D:
    sz = size(slice);
    if isvector(slice) || numel(sz)==2 && (sz(1)==1 || sz(2)==1)
        sliceType = 'vector';
        y = mapComplex(slice, dispOpt);
        h = plot(ax, y, 'LineWidth',1.5);
        xlabel(ax, 'Index');
        ylabel(ax, getYLabel(dispOpt));
        grid(ax,'on');
        hPlot = h;
    elseif ismatrix(slice)  % 2D: sz(1)>1 und sz(2)>1
        sliceType = '2D';
        img = mapComplex(slice, dispOpt);
        hImg = imagesc(ax, img);
        axis(ax,'image');
        if ~isempty(cmap)
            colormap(ax, cmap);
        end
        if isempty(clim)
            set(ax, 'CLimMode','auto');
        else
            caxis(ax, clim);
        end
        colorbar(ax);
        xlabel(ax, '');
        ylabel(ax, '');
        hPlot = hImg;
    else
        error(['Slice ist mehrdimensional (>2D) nach squeeze. ', ...
               'Bitte vorher in 2D reduzieren oder nur 1D/2D-Slices verwenden.']);
    end
end

% Hilfsfunktion: Wandle komplexe slice-Werte gemäß dispOpt um
function out = mapComplex(slice, dispOpt)
    if isreal(slice)
        out = slice;
        return;
    end
    switch dispOpt
        case 'abs'
            out = abs(slice);
        case 'real'
            out = real(slice);
        case 'imag'
            out = imag(slice);
        case 'angle'
            out = angle(slice);
        otherwise
            % Default: Bei 2D-Plot: Betrag; bei 1D-Plot: Realteil
            if ismatrix(slice) && size(slice,1)>1 && size(slice,2)>1
                out = abs(slice);
            else
                out = real(slice);
            end
    end
end

% Hilfsfunktion: Y-Achsen-Beschriftung je nach dispOpt
function lbl = getYLabel(dispOpt)
    switch dispOpt
        case 'abs'
            lbl = '|Amplitude|';
        case 'real'
            lbl = 'Real';
        case 'imag'
            lbl = 'Imag';
        case 'angle'
            lbl = 'Phase (rad)';
        otherwise
            lbl = 'Amplitude';
    end
end
