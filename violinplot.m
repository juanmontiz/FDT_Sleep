function violins = violinplot(data, cats, varargin)
%Violinplots plots violin plots of some data and categories
%   VIOLINPLOT(DATA) plots a violin of a double vector DATA
%
%   VIOLINPLOT(DATAMATRIX) plots violins for each column in
%   DATAMATRIX.
%
%   VIOLINPLOT(DATAMATRIX, CATEGORYNAMES) plots violins for each
%   column in DATAMATRIX and labels them according to the names in the
%   cell-of-strings CATEGORYNAMES.
%
%   In the cases above DATA and DATAMATRIX can be a vector or a matrix,
%   respectively, either as is or wrapped in a cell.
%   To produce violins which have one distribution on one half and another
%   one on the other half, DATA and DATAMATRIX have to be cell arrays
%   with two elements, each containing a vector or a matrix. The number of
%   columns of the two data sets has to be the same.
%
%   VIOLINPLOT(DATA, CATEGORIES) where double vector DATA and vector
%   CATEGORIES are of equal length; plots violins for each category in
%   DATA.
%
%   VIOLINPLOT(TABLE), VIOLINPLOT(STRUCT), VIOLINPLOT(DATASET)
%   plots violins for each column in TABLE, each field in STRUCT, and
%   each variable in DATASET. The violins are labeled according to
%   the table/dataset variable name or the struct field name.
%
%   violins = VIOLINPLOT(...) returns an object array of
%   <a href="matlab:help('Violin')">Violin</a> objects.
%
%   VIOLINPLOT(..., 'PARAM1', val1, 'PARAM2', val2, ...)
%   specifies optional name/value pairs for all violins:
%     'Width'        Width of the violin in axis space.
%                    Defaults to 0.3
%     'Bandwidth'    Bandwidth of the kernel density estimate.
%                    Should be between 10% and 40% of the data range.
%     'ViolinColor'  Fill color of the violin area and data points. Accepts
%                    1x3 color vector or nx3 color vector where n = num
%                    groups. In case of two data sets being compared it can 
%                    be an array of up to two cells containing nx3
%                    matrices.
%                    Defaults to the next default color cycle.
%     'ViolinAlpha'  Transparency of the violin area and data points.
%                    Can be either a single scalar value or an array of
%                    up to two cells containing scalar values.
%                    Defaults to 0.3.
%     'MarkerSize'   Size of the data points, if shown.
%                    Defaults to 24
% 'MedianMarkerSize' Size of the median indicator, if shown.
%                    Defaults to 36
%     'EdgeColor'    Color of the violin area outline.
%                    Defaults to [0.5 0.5 0.5]
%     'BoxColor'     Color of the box, whiskers, and the outlines of
%                    the median point and the notch indicators.
%                    Defaults to [0.5 0.5 0.5]
%     'MedianColor'  Fill color of the median and notch indicators.
%                    Defaults to [1 1 1]
%     'ShowData'     Whether to show data points.
%                    Defaults to true
%     'ShowNotches'  Whether to show notch indicators.
%                    Defaults to false
%     'ShowMean'     Whether to show mean indicator
%                    Defaults to false
%     'ShowBox'      Whether to show the box.
%                    Defaults to true
%     'ShowMedian'   Whether to show the median indicator.
%                    Defaults to true
%     'ShowWhiskers' Whether to show the whiskers
%                    Defaults to true
%     'GroupOrder'   Cell of category names in order to be plotted.
%                    Defaults to alphabetical ordering

% Copyright (c) 2016, Bastian Bechtold
% This code is released under the terms of the BSD 3-clause license

% THIRD-PARTY FILE — not included in this repository.
% Download from: https://github.com/bastibe/Violinplot-Matlab
error('violinplot.m is a third-party file. Download from: https://github.com/bastibe/Violinplot-Matlab');
end
