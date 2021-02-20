%-------------------------------------------------------------------------
% LTF Method for leaflet Ka moduli
%
% Method based on Doktorova et al. 2019, Biophysical Journal 116, 487-502
%
% Milka Doktorova, March 2019
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% get_bilayerKa_stats.m
%
% Calculates bilayer Ka as the harmonic average of the two leaflet Ka-s
% If leaflet Ka errors are provided, calculates and returns the propagated
% error on the bilayer Ka as well
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Input parameters
%
% topKa     Ka of top leaflet
% botKa     Ka of bottom leaflet
% topKaerr  error of Ka of top leaflet (optional)
% botKaerr  error of Ka of bottom leaflet (optional)
%
% Output
% 
% bilayer-Ka OR [bilayer-Ka bilayer-Ka-error]
%-------------------------------------------------------------------------
function [out] = get_bilayerKa_stats(topKa, botKa, topKaerr, botKaerr)

bKa = 1/(0.5*(1/topKa+1/botKa));

if ~exist('topKaerr','var') || ~exist('botKaerr','var')
    out = bKa;
else 
    bKaerr = 2*sqrt((topKaerr^2)/(topKa^4) + (botKaerr^2)/(botKa^4))/((1/topKa+1/botKa)^2);
    out = [bKa bKaerr];
end

end