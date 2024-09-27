%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gerchberg-Saxton Algorithm for Phase Retrieval and SLM Projection
% Author: Hisay Lama
% Email: hisaylama@gmail.com
% 
% This code uses the Gerchberg-Saxton algorithm to retrieve the phase of
% a desired intensity pattern and project it onto an SLM. The projection 
% and display functions are borrowed from 'ots1.0.1' by Volpe et al.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace, close all figures
clear; clc; close all;

%% Input Parameters
M = 800; % SLM side length (pixels)
N = 600; % SLM side height (pixels)
psize = 20e-6; % Pixel size (meters)
lambda = 1064e-9; % Laser wavelength (meters)

%% Coordinate System
x = ([0.5:1:M-0.5] - M/2) * psize;
y = ([0.5:1:N-0.5] - N/2) * psize;
[X, Y] = meshgrid(x, y);

%% Input Beam Parameters
w_beam = 250 * psize; % Radius of the input beam (pixels)
u_in = circ(sqrt(X.^2 + Y.^2) / w_beam); % Input beam amplitude

%% Display Input Beam Irradiance
I_in = abs(u_in.^2); % Intensity of input beam
figure(1), imagesc(x, y, I_in), colormap('gray');
xlabel('x (m)'), ylabel('y (m)');
axis xy;

%% Generate Random Phase
random_phase = 2 * pi * rand([N, M]); 

%% Get Target Pattern from User
num_spots = input('Enter the number of spots: ');
spot_separation = input('Enter the inter-spot separation (pixels): ');
u_target = generate_spot_pattern(M, N, num_spots, spot_separation * psize);

%% Display Target Pattern
I_target = abs(u_target.^2);
figure(2), imagesc(x, y, I_target), colormap('gray');
xlabel('x (m)'), ylabel('y (m)');
axis xy;

%% GS Algorithm Parameters
focal_length = input('Enter focal length (pixels): '); % Focal length in pixels
z_focus = 0; % Distance from focus
num_iterations = input('Maximum number of iterations: ');

%% GS Algorithm - Phase Retrieval
k = 2 * pi / lambda; % Wave number
[fx, fy] = meshgrid((1 / psize) * [-M/2:1:M/2-1] / M, (1 / psize) * [-N/2:1:N/2-1] / N);
x1 = lambda * focal_length * fx;
y1 = lambda * focal_length * fy;

% Start iterations to retrieve phase mask
Phase_hot = random_phase; % Initial random phase
for iter = 1:num_iterations
    Efocus = exp(1i * k * (2 * focal_length + z_focus + x1.^2 + y1.^2)) / (1i * lambda * focal_length) ...
             .* fftshift(fft2(abs(u_in) .* exp(-1j * Phase_hot) .* exp((-1i * pi * z_focus / (lambda * focal_length^2)) * (X.^2 + Y.^2))));
    phase_focus = mod(angle(Efocus), 2 * pi); 
    E_target = abs(u_target) .* exp(-1j * phase_focus);
    E_inverse = ifftshift(ifft2(ifftshift(E_target)));
    Phase_hot = mod(angle(E_inverse), 2 * pi);
end

%% Generate Phase Mask for SLM
final_phase_mask = mod(Phase_hot, 2 * pi);
figure(3), imshow(mat2gray(final_phase_mask + 1.2 * pi)), axis equal tight off;

%% Projecting Phase Mask onto SLM
project_on_slm(final_phase_mask, M, N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u_target = generate_spot_pattern(M, N, num_spots, separation)
    % Generates a target pattern with multiple spots.
    % Inputs:
    %   M - SLM width in pixels
    %   N - SLM height in pixels
    %   num_spots - number of spots
    %   separation - distance between spots in meters
    %
    % Outputs:
    %   u_target - target amplitude pattern
    
    % Define the spot radius
    spot_radius = 1 * 20e-6; % Example spot size
    
    % Initialize target pattern
    u_target = zeros(N, M);
    
    % Place spots based on user input
    for i = 1:num_spots
        % Place the spot at a defined separation
        x_center = M/2 + (i-1) * separation;
        y_center = N/2 + (i-1) * separation;
        
        % Apply a circular spot with the defined radius
        u_target = u_target + circ(sqrt((x_center - (1:M)).^2 + (y_center - (1:N)).^2) / spot_radius);
    end
end

function project_on_slm(phase_mask, M, N)
    % Projects the phase mask onto a Spatial Light Modulator (SLM).
    % Inputs:
    %   phase_mask - calculated phase mask
    %   M - SLM width in pixels
    %   N - SLM height in pixels
    
    % Determine monitor setup
    scrnum = 2; % Assuming the SLM is the second monitor
    ScreenPos = get(0, 'MonitorPositions');
    ScreenWidth = ScreenPos(scrnum, 3) - ScreenPos(scrnum, 1) + 1;
    ScreenHeight = ScreenPos(scrnum, 4) - ScreenPos(scrnum, 2) + 1;
    
    % Ensure screen size matches SLM size
    if [M, N] ~= [ScreenHeight, ScreenWidth]
        warning('The size of the phase mask and screen must match!');
    end
    
    % Create figure for SLM display
    slmfig = figure('Name', 'SLM Screen', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                    'ToolBar', 'none', 'Units', 'pixels', 'Resize', 'off', ...
                    'Position', [ScreenPos(scrnum, 1), ScreenPos(scrnum, 2), ScreenWidth, ScreenHeight]);
    axes('Units', 'Normalized', 'Position', [0, 0, 1, 1], 'Visible', 'off');
    image(phase_mask / (2 * pi) * 255);
    axis equal tight off;
    colormap(gray(255));
    drawnow();
    
    % Close other figures, keep only SLM display
    fh = findall(0, 'type', 'figure').';
    close(fh(fh ~= slmfig));
end
