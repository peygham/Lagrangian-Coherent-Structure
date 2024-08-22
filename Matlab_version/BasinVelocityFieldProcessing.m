%% Code for Velocity Data Processing
% this code reads FVCOM results and make a gridded velocity field for LCS
% you need M.mat file of FVCOM, otherwise adjust the code to read lat/lon
% from the output nc file

% phg@akvaplan.niva.no

% Clear workspace, close all figures, and clear command window
close all; clear; clc;

% Load the mesh data from the specified path
load('F:\RFF_NordLand\Model\Nland_v3\mat_files\M.mat');

% Define the path where the NetCDF files are stored
filepath = 'Y:\nird\projects\NS9067K\apn_backup\FVCOM\FVCOM_results_hansi\nordL_2023.Nordland_2023\output03_2023\';
filepattern = fullfile(filepath, 'Nordland_2023_00*.nc');
thefiles = dir(filepattern); % List all matching files
file_no = length(thefiles); % Count the number of files

%% Extract and Combine Time and Date Information
mydate = []; % Initialize an empty array to store time and date

% Loop through the specified file indices to extract time and date
for ii = 23:23 % You can modify this range based on your needs
    filepathname = fullfile(filepath, thefiles(ii).name); % Full path to the file
    % Read and convert time information from the NetCDF file
    tem_Itime = double(ncread(filepathname, 'Itime')) + datenum(1858, 11, 17, 0, 0, 0);
    tem_Itime2 = double(ncread(filepathname, 'Itime2'));
    tem_date = tem_Itime + tem_Itime2 / (3600000 * 24); % Convert to MATLAB date format
    mydate = [mydate; tem_date]; % Concatenate the dates
end
clearvars tem_Itime tem_Itime2 tem_date; % Clear temporary variables

%% Define the Basin Boundary
% Convert UTM coordinates to latitude and longitude for mesh elements
[Lat, Lon] = utm2deg(Mobj.xc, Mobj.yc, '33 W');

% Load basin definition data (assuming it contains lonv, latv, and related variables)
load('K.mat'); % Load the K.mat file containing basin and island data

% Determine which mesh elements fall within the defined basin boundary
in = inpolygon(Lon, Lat, K.lonv, K.latv);
inBasin = find(in > 0); % Indices of elements within the basin
mylon = Lon(inBasin); % Longitudes of elements within the basin
mylat = Lat(inBasin); % Latitudes of elements within the basin

% Create a combined array of basin coordinates
basin = [mylon, mylat];
basinElems = find(in > 0); % Indices of basin elements

% Load island boundary definitions
load('IsK.mat'); % Load IsK.mat file containing island boundary data

% Extract basin boundary coordinates from the basin definition
lonB = basin(K.k, 1);
latB = basin(K.k, 2);

%% Grid Setup for Meshing
% Define the resolution of the grid (can be adjusted based on computational resources)
grid_res = 5000;
myxm = linspace(min(K.lonv), max(K.lonv), grid_res);
myym = linspace(min(K.latv), max(K.latv), grid_res);
[X, Y] = meshgrid(myxm, myym); % Create a grid for interpolation

% Determine which grid points fall inside the basin and islands
inB = inpolygon(X, Y, lonB, latB);
inIs1 = inpolygon(X, Y, IsK.lonIs1, IsK.latIs1);
inIs2 = inpolygon(X, Y, IsK.lonIs2, IsK.latIs2);
inIs3 = inpolygon(X, Y, IsK.lonIs3, IsK.latIs3);
inIs4 = inpolygon(X, Y, IsK.lonIs4, IsK.latIs4);
inIs5 = inpolygon(X, Y, IsK.lonIs5, IsK.latIs5);
inIs6 = inpolygon(X, Y, IsK.lonIs6, IsK.latIs6);

%% Velocity Data Processing
% Loop through each NetCDF file to process velocity data
for ii = 1:file_no
    disp(['Processing file: ', num2str(ii)]); % Display progress
    
    u = []; % Initialize u velocity array
    v = []; % Initialize v velocity array
    filepathname = fullfile(filepath, thefiles(ii).name); % Full path to the file
    
    % Loop through each hourly timestep in the file
    for jj = 1:24
        disp(['Processing hour: ', num2str(jj)]); % Display progress
        
        % Read u and v velocity components for the current hour
        u_tem = double(ncread(filepathname, 'u', [1 1 jj], [Inf 1 1]));
        v_tem = double(ncread(filepathname, 'v', [1 1 jj], [Inf 1 1]));
        
        % Extract velocities corresponding to the basin elements
        u_basin = u_tem(inBasin);
        v_basin = v_tem(inBasin);
        
        % Concatenate velocity data across all hours
        if jj == 1
            u = u_basin;
            v = v_basin;
        else
            u = [u, u_basin];
            v = [v, v_basin];
        end
    end
    
    % Meshing and Interpolating Data
    U = zeros(size(u, 2), size(X, 1), size(X, 2)); % Pre-allocate U array
    V = zeros(size(v, 2), size(X, 1), size(X, 2)); % Pre-allocate V array
    
    % Loop through each time step to interpolate the velocity field
    for kk = 1:size(u, 2)
        % Interpolate u and v velocities onto the grid
        Fu = scatteredInterpolant(mylon, mylat, double(u(:, kk)));
        Fv = scatteredInterpolant(mylon, mylat, double(v(:, kk)));
        uq = Fu(X, Y);
        uq(~inB | inIs1 | inIs2 | inIs3 | inIs4 | inIs5 | inIs6) = nan; % Set values outside basin and islands to NaN
        vq = Fv(X, Y);
        vq(~inB | inIs1 | inIs2 | inIs3 | inIs4 | inIs5 | inIs6) = nan; % Set values outside basin and islands to NaN
        U(kk, :, :) = uq; % Store interpolated u velocity
        V(kk, :, :) = vq; % Store interpolated v velocity
    end
    
    % Save interpolated velocity fields to a .mat file
    field.u = U;
    field.v = V;
    save(fullfile(filepath, ['VelocityMar5_' num2str(ii, '%02d') '.mat']), 'field');
end

%% Save Final Mesh and Time Data
S.X = X; % Save X grid coordinates
S.Y = Y; % Save Y grid coordinates
S.time = mydate; % Save time data
save(fullfile(filepath, 'MarMesh5.mat'), 'S', '-v7.3'); % Save the mesh and time data to a .mat file
