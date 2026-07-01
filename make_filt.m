function [bfilt, afilt] = make_filt(flp,fhi,delt,k)
% MAKE_FILT Designs a Butterworth bandpass filter for use with filtfilt.
% Author: Juan Manuel Monti
%
% Inputs:
%   flp   lower cutoff frequency (Hz), e.g. 0.008
%   fhi   upper cutoff frequency (Hz), e.g. 0.08
%   delt  sampling interval (s), typically TR
%   k     filter order (e.g. 2 for 2nd-order Butterworth)
%
% Outputs:
%   bfilt  numerator coefficients for filtfilt
%   afilt  denominator coefficients for filtfilt

fnq = 1/(2*delt);       % Nyquist frequency
Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency

[bfilt, afilt] = butter(k,Wn);   % construct the filter
end
