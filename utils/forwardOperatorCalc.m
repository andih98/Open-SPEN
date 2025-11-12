function [kernel, ky, s] = forwardOperatorCalc(kspace, param,use)
% FORWARDOPERATORCALC Calculates the analytical SPEN forward operator kernel.
%
% This function determines the quadratic phase modulation kernel for SPEN,
% which relates object space (y) to frequency space (ky) over time (s).
%
% Args:
%   kspace (struct): Struct containing k-space trajectory data (e.g., ktraj_adc).
%   param (struct): Struct containing sequence parameters (R, fov, Ny, etc.).
%
% Returns:
%   kernel (matrix): The calculated SPEN forward operator kernel (time x space).
%   ky (vector): The k-space values used in the calculation.
%   s (vector): The effective spatial shifts corresponding to the k-space values.

    R = param.R;        % Acceleration factor
    fov = param.fov(1); % Field of View in the SPEN dimension [m]
    Ny = param.Ny;      % Number of phase encoding steps

    % Calculate the constant 'a' for the quadratic phase term: exp(i * a * (y - s)^2)
    % The original logic for 'a' calculation is complex but preserved here.
    % a = 2 * R * G^2 / (2 * (G * FOV / T_rf))
    % This simplifies to a constant incorporating system parameters and sequence timings.
    a = 2 * (param.sweepBw * param.rfdur / param.Ny) * param.g^2 / (2 * (param.g * fov / param.rfdur));

    %% --- Determine k-space sampling points (ky) ---

    % Logic for sequences with a SPEN navigator
    ky_adc = kspace.ktraj_adc(2, :);
    % Find indices where k-space trajectory changes (start of new acquisition lines)
    idx_change = find(diff(ky_adc) ~= 0);
    idx_start = [1, idx_change + 1];
    ky_lines = ky_adc(idx_start); % Extract ky values at the start of each line

    if param.spenNavigator
     
        if size(ky_lines, 1) < 2
            % If only one row (common for some k-space calculations), flip for correct order
            ky = fliplr(ky_lines);
        else
            % Sort the ky values to get the correct order for the SPEN dimension
            ky = sort(ky_lines(2:Ny + 1));
        end
    else
        % Logic for standard SPEN (no explicit navigator block)
        ky = sort(ky_lines(1:Ny));
    end
    %% --- Define Spatial and Shift Vectors ---

    % Spatial vector 'y' covering the FOV
    % The number of spatial points is set by the time-Bandwidth product (T_rf * BW_sweep)
     if strcmp(use,'exc')
        y = linspace(-fov / 2, fov / 2, round(abs(param.sweepBw * param.rfdur)));
    elseif strcmp(use,'ref')
        y = linspace(-fov , fov , 2*round(abs(param.sweepBw * param.rfdur)));
    else
        warning('Define SPEN mode "exc" when the chirped RF is used as excitation pulse or "ref" when it is used for refocusing!')
        return
    end
    % Effective spatial shift vector 's' corresponding to k-space points (ky)
    % This shift s = (G * T_k * FOV) / (G_max * T_rf) is proportional to ky
    s = fov * ky / (R * Ny / fov);

    %% --- Calculate Quadratic Phase Kernel ---

    % The kernel is built by iterating over each time/ky point (s(i))
    % The phase is calculated as Phi = a * (y - s)^2
    kernel = zeros(size(s, 2), size(y, 2)); % Preallocate matrix (Time/ky x Space)
    for i = 1:size(s, 2)
        % Calculate the phase term for the current k-space point/shift 's(i)'
        kernel(i, :) = a * (y - s(i)).^2;
    end

    %% --- Final Kernel Calculation ---

    % 1. Apply the exponential to get the complex phase modulation: exp(i * Phi)
    % 2. Apply FFT along the spatial dimension (axis 2) to transform the spatial phase
    %    into a frequency-dependent convolution kernel.
    kernel = ifftshift(fft(fftshift(exp(-1i * kernel), 2), [], 2), 2);

end

