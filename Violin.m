classdef Violin < handle
    % Violin creates violin plots for some data
    %   A violin plot is an easy to read substitute for a box plot
    %   that replaces the box shape with a kernel density estimate of
    %   the data, and optionally overlays the data points itself.
    %   It is also possible to provide two sets of data which are supposed
    %   to be compared by plotting each column of the two datasets together
    %   on each side of the violin.
    %
    %   Additional constructor parameters include the width of the
    %   plot, the bandwidth of the kernel density estimation, the
    %   X-axis position of the violin plot, and the categories.
    %
    %   Use <a href="matlab:help('violinplot')">violinplot</a> for a
    %   <a href="matlab:help('boxplot')">boxplot</a>-like wrapper for
    %   interactive plotting.
    %
    %   See for more information on Violin Plots:
    %   J. L. Hintze and R. D. Nelson, "Violin plots: a box
    %   plot-density trace synergism," The American Statistician, vol.
    %   52, no. 2, pp. 181-184, 1998.
    %
    % Violin Properties:
    %    ViolinColor    - Fill color of the violin area and data points.
    %                     Can be either a matrix nx3 or an array of up to two
    %                     cells containing nx3 matrices.
    %                     Defaults to the next default color cycle.
    %    ViolinAlpha    - Transparency of the violin area and data points.
    %                     Can be either a single scalar value or an array of
    %                     up to two cells containing scalar values.
    %                     Defaults to 0.3.
    %    EdgeColor      - Color of the violin area outline.
    %                     Defaults to [0.5 0.5 0.5]
    %    BoxColor       - Color of the box, whiskers, and the outlines of
    %                     the median point and the notch indicators.
    %                     Defaults to [0.5 0.5 0.5]
    %    MedianColor    - Fill color of the median and notch indicators.
    %                     Defaults to [1 1 1]
    %    ShowData       - Whether to show data points.
    %                     Defaults to true
    %    ShowNotches    - Whether to show notch indicators.
    %                     Defaults to false
    %    ShowMean       - Whether to show mean indicator.
    %                     Defaults to false
    %    ShowBox        - Whether to show the box.
    %                     Defaults to true
    %    ShowMedian     - Whether to show the median indicator.
    %                     Defaults to true
    %    ShowWhiskers   - Whether to show the whiskers
    %                     Defaults to true
    %    HalfViolin     - Whether to do a half violin(left, right side) or
    %                     full. Defaults to full.
    %    QuartileStyle - Option on how to display quartiles, with a
    %                     boxplot, shadow or none. Defaults to boxplot.
    %    DataStyle      - Defines the style to show the data points. Opts: 
    %                     'scatter', 'histogram' or 'none'. Default is 'scatter'.
    %
    % Violin Children:
    %    ScatterPlot    - <a href="matlab:help('scatter')">scatter</a> plot of the data points
    %    ScatterPlot2   - <a href="matlab:help('scatter')">scatter</a> second plot of the data points
    %    ViolinPlot     - <a href="matlab:help('fill')">fill</a> plot of the kernel density estimate
    %    ViolinPlot2    - <a href="matlab:help('fill')">fill</a> second plot of the kernel density estimate
    %    BoxPlot        - <a href="matlab:help('fill')">fill</a> plot of the box between the quartiles
    %    WhiskerPlot    - line <a href="matlab:help('plot')">plot</a> between the whisker ends
    %    MedianPlot     - <a href="matlab:help('scatter')">scatter</a> plot of the median (one point)
    %    NotchPlots     - <a href="matlab:help('scatter')">scatter</a> plots for the notch indicators
    %    MeanPlot       - line <a href="matlab:help('plot')">plot</a> at mean value

    % Copyright (c) 2016, Bastian Bechtold
    % This code is released under the terms of the BSD 3-clause license

    % THIRD-PARTY FILE — not included in this repository.
    % Download from: https://github.com/bastibe/Violinplot-Matlab

    methods
        function obj = Violin(varargin)
            error('Violin.m is a third-party file. Download from: https://github.com/bastibe/Violinplot-Matlab');
        end
    end
end
