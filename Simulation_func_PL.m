% this file generates the average of 1000 point loss variable instances for
% a given propagation distance, laser wavelength, transmit beamwidth and
% receive aperture

% Written by: Mateen Ashraf

function [gains] = Simulation_func_PL(Prop_dist, r_aper, w_0, lambda_light)

epsilon = 6.7; % 6.7 corresponds to negligible pointing error, 1.1 corresponds to severe pointing error

N  = 1000; % number of realizations

variable_gain = zeros(1,N); % dummy variable

for i = 1:N

    variable_gain(i) = pointing_error_rv(Prop_dist*10^3,w_0,lambda_light,epsilon,r_aper); % propagation distance is assumed to be in kilo meters here

end

gains = mean(variable_gain); % just the pointing loss average over multiple instances/realizations of the pointing loss variable 