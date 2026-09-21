% This function generates one instance of point loss for a given
% propagation distance, wavelength, visibility factor and receive aperture.

% Written by: Mateen Ashraf

% These formulas can be found in 
% [1] T. V. Nguyen, H. D. Le, N. T. Dang and A. T. Pham, "Average Transmission Rate and Outage Performance of Relay-Assisted Satellite Hybrid FSO/RF Systems," 2021 International Conference on Advanced Technologies for Communications (ATC), Ho Chi Minh City, Vietnam, 2021, pp. 1-6.
% [2] F. Yang, J. Cheng and T. A. Tsiftsis, "Free-Space Optical Communication with Nonzero Boresight Pointing Errors," in IEEE Transactions on Communications, vol. 62, no. 2, pp. 713-725, February 2014.

function x = pointing_error_rv(z, w0, lambda, epsilon, r)

    wz = w0 * sqrt(1 + (lambda * z / (pi * w0^2))^2); % see the references for these formulas

    A0 = (erf( sqrt(pi)*r / (sqrt(2)*wz) ))^2; % see the references for these formulas

    u = rand(); % uniform random variable generation

    x = A0 * u^(1/(epsilon^2)); % power law distribution

end
