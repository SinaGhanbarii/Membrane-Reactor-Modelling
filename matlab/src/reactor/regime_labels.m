function code = regime_labels(beta, pinch, flux_on, idle)
% Local regime code along the bed (a proposal; thresholds to be agreed):
%   1 no flux     retentate H2 partial pressure is not above the permeate pressure
%   2 reaction    beta < 0.9, the reaction is not at equilibrium
%   3 pinch       outlet H2 partial pressure within 5 % of the permeate pressure
%   4 permeance   reaction at equilibrium and driving force left, so permeation is limited by the membrane
%   5 idle        reaction rate and H2 flux are both below 1 % of their maximum along the bed: nothing is limiting
code = 4 * ones(size(beta));
code(pinch < 1.05) = 3;
code(beta < 0.9) = 2;
code(~flux_on) = 1;
if nargin > 3, code(idle) = 5; end
end
