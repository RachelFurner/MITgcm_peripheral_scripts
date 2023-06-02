%% ------------------------------------------------------------------------
% DRM, 04/08/11.
% m-file: gen_so_inputs.m

% This m-file is for generating the input files for a channel model of the
% ACC including a ridge and two (optional) obstacles, one on the northern
% boundary and one on the southern boundary. Restoring zones to the north and
% south can also be separatley activated and deactivated.

%% ------------------------------------------------------------------------
% Flush the memory.

clc;
clear;
close all;

%% ------------------------------------------------------------------------
% Make some specifications of the problem Im going to create my inputs for.
% Breakdown of the file path into units.
name.directory = './';
name.problem = 'so-010km-jet/flatb/';

% Choose whether to output the inputs or not.
output = 1;

% Choose whether files should be output as big-endian or little-endian.
format = 'b';

% Select bathymetry options.
ridge = 0;
north = 0;
south = 0;

%% ------------------------------------------------------------------------
% Problem dependent parameters.

% Choose the amplitude of the temperature variations.
dtheta = 7.;

% Choose the wind stress magnitude.
tau0 = 0.20;
tau1 = 0.0;

% Choose the amplitude of the surface heat flux.
hflux0 = 10.;

% Choose the e-folding scale for the stratification.
z0 = 1000.;

% Half-width of the ridge.
Lr = 500.E3;

% Height of the ridge.
Hr = 0.;

% Latitudinal extent of the barriers.
LbN = 0.   %LbN = 1000.E3; % Barrier extending south from Northern boundary.
LbS = 0.   %LbS = 1000.E3; % Barrier externing north from Southern boundary,

%% ------------------------------------------------------------------------
% Grid spacing in degrees and the number of vertical levels.
dx = 10.E3;
dy = dx;
dz = [ 10.0000, 11.3752, 12.9394, 14.7188, 16.7429, 19.0453, 21.6643, ...
    24.6435, 28.0323, 31.8872, 36.2722, 41.2602, 46.9342, 53.3883, ...
    60.7301, 69.0814, 78.5812, 89.3874, 101.6795, 115.6621, 131.5674, ...
    149.6600, 170.2406, 193.6514, 220.8552, ...
    250.0000, 250.0000, 250.0000, 250.0000, 250.0000, ...
    250.0000, 250.0000, 250.0000, 250.0000, ...
    250.0000, 250.0000, 250.0000, 250.0000 ];
  
% Specify the number of extra grid points to have outside of the channel.
kx = 0;
ky_north = 3;
ky_south = 3;

% pi to machine precision.
pi = 4.0*atan( 1.0 );

%% ------------------------------------------------------------------------
% Specify the parameters that define the channel dimensions.

% Length of the channel.
%Lx = 9600.E3;
Lx = 2400.E3;

% Width of the channel.
%Ly = 2000.E3;
Ly = 980.E3;

% Depth of the channel.
Lz = 5000.;

% Set the number of grid boxes.
nx = Lx/dx + kx;
ny = Ly/dy + ky_north + ky_south;
nz = 38;

% Calculate the x and y locations of the tracer points.
x = repmat( dx*( (1:1:nx)'-0.5 ) - 0.5*Lx, [ 1 ny nz ] );
y = repmat( dy*( (1:1:ny)-0.5 ) - 0.5*Ly-ky_south*dy, [ nx 1 nz ] );

% Calculate the Coriolis frequency, normalised by f0.
f = ( -1.11E-4 + 1.47E-11*y(:,:,1) ) / -1.11E-4;

% Calculate the depth of the tracer points.
z = NaN( nx, ny, nz );

% Set depth of first level.
z( :, :, 1 ) = 0.5*dz( 1 );

% Loop over remaining depths.
for k = 2:1:nz;
    z( :, :, k ) = z( :, :, k-1 ) + 0.5*( dz( k-1 ) + dz( k ) );
end;

% Set the depth to negative.
z = -z;

%% ------------------------------------------------------------------------
% Define the bathymetry.

% Begin with a flat bottom at 5000 m.
depth = -Lz*ones( nx, ny );

% Use the bathymetry flags to set any obstacles.

% Add ridge across centre of channel.
if ridge    
    for j = 1:1:ny;
        for i = 1:1:nx
            if x( i, j, 1 ) > -Lr && x( i, j, 1 ) < Lr
                depth( i, j ) = -Lz ...
                    + 0.5*Hr*( 1. + cos( pi*x( i, j, 1 )/Lr ) );
            else
                depth( i, j ) = -Lz;
            end;
        end;
    end;
end;

% Add continent to southern boundary.
if south
    by = LbS/dy;
    depth( 231:250, 1:by ) = 0.;
end;

% Add continent to the northern boundary.
if north
    by = ( Ly - LbN )/dy;
    depth( 711:730, by:end ) = 0.; % This line gives the northern boundary ridge 1 grid box more than it should be in y, but is consistent with so-010km-diffkr expts.
%    depth( 711:730, by+1:end ) = 0.; % This line makes the northern boundary ridge the correct extent, but is inconsistent with so-010km-diffkr expts.
end;

% Put a solid wall on the northern boundary.
depth( :, end-ky_north+1:end ) = 0.;

% Put a solid wall on the southern boundary.
depth( :, 1:ky_south ) = 0.;
%% ------------------------------------------------------------------------
% Generate the spatially varying restoring mask.

mask = zeros(nx,ny,nz);

% Set the surface layer to restoring everywhere.
mask( :, :, 1 ) = 1.;

%% ------------------------------------------------------------------------
% Generate the idealised temperature stratification, used as both an initial
% condition and for the restoring temperature in the sponge regions.

% Linear gradient at surface with exponential decay at depth and 0oC at
% the southern boundary (similar to Abernathey et al., 2011).
t = dtheta * ( (y+.5*Ly)/Ly ) .* ( exp( z/z0 ) - exp( -Lz/z0 ) ) / ...
    ( 1. - exp( -Lz/z0 ) );

%% ------------------------------------------------------------------------
% Construct an initial temperature distribution with a little bit of white
% gaussian noise added; this helps encourage instability.

% Preallocate for speed.
noisyt = NaN( size( t ) );

% Add some white Gaussian noise to the initial temperature.
for k = 1:1:nz;
    noisyt( :, :, k ) = awgn( t( :, :, k ), 30 );
end;

% Make sure there is no water colder than 0.5oC.
%noisyt( noisyt < 0.25 ) = 0.25;

%% ------------------------------------------------------------------------
% Calculate the idealised wind stress profile.

tau = tau0*0.5*( 1. + cos( 2.0*pi*y( :, :, 1 )/Ly ) ) + tau1;
%tau = f.*tau0.*ones( nx, ny );

tauy = 0.25*tau0*exp( -0.5*x( :, :, 1 ).*x( :, :, 1 )/250000.E6 );

%% ------------------------------------------------------------------------
% Generate surface heat flux.

hflux = -hflux0*sin( pi*y( :, :, 1 )/Ly );

% ------------------------------------------------------------------------
% Generate the spatially varying diapycnal diffusivity.

diffkr = 1.E-5*ones(nx,ny,nz);

% Increase the diapycnal diffusivity in the northern sponge regions.
%peak_diff = 5.E-3
peak_diff = 10.E-3
%increased_diff_length = 150.E3
increased_diff_length = 75.E3

disp(Ly/2 - increased_diff_length)

for j = 1:1:ny;
  if y( 1, j ) > (Ly/2 - increased_diff_length)  % Want peak at land edge (490)
      %disp(y(1,j))
      diffkr( :, j, : ) = 1.E-5 ...
          + 0.5*peak_diff*( 1 + cos( pi*( y(:,j,:) - .5*Ly )/increased_diff_length ) ) ...
          - 0.5*1.E-5*( 1 + cos( pi*( y(:,j,:) - .5*Ly )/increased_diff_length ) );
      %disp(diffkr(1,j,1))
  end;
end;

%% ------------------------------------------------------------------------

% Draw a map of the bathymetry.
figure;
pcolor( x( :, :, 1 ), y ( :, :, 1 ), depth );
shading flat;
axis equal tight;
colorbar;
title( 'Bathymetry' );

%% ------------------------------------------------------------------------

% Draw a map of the surface forcing.
figure;

subplot( 211 );
pcolor( x( :, :, 1 ), y ( :, :, 1 ), t( :, :, 1 ) );
shading flat;
axis equal tight;
colorbar;
title( 'Surface restoring \theta' );

subplot( 212 );
pcolor( x( :, :, 1 ), y ( :, :, 1 ), noisyt( :, :, 1 )-t( :, :, 1 ) );
shading flat;
axis equal tight;
colorbar;
title( 'Surface initial \theta' );

figure;

subplot( 211 );
pcolor( x( :, :, 1 ), y ( :, :, 1 ), tau( :, :, 1 ) );
shading flat;
axis equal tight;
colorbar;
title( '\tau' );

subplot( 212 );
pcolor( x( :, :, 1 ), y ( :, :, 1 ), diffkr( :, :, 1 ) );
shading flat;
axis equal tight;
colorbar;
title( 'diffkr' );

%% ------------------------------------------------------------------------
% Read out the depths into a binary bathymetry file.

if output;
    %out_file = fullfile( name.directory, name.problem, 'input/bathy.bin' );
    out_file = 'bathy.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( depth, nx*ny, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;

%% ------------------------------------------------------------------------
% Read out the restoring mask to binary file.

if output;
    %out_file = fullfile( name.directory, name.problem, 'input/mask.bin' );
    out_file = 'mask.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( mask, nx*ny*nz, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;

%% ------------------------------------------------------------------------
% Read out the diapycnal diffusivity map to binary file.

if output;
    %out_file = fullfile( name.directory, name.problem, 'input/diffkr.bin' );
    out_file = 'diffkr.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( diffkr, nx*ny*nz, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;

%% ------------------------------------------------------------------------
% Read out the restoring temperature to binary file.

if output;
    %out_file = fullfile( name.directory, name.problem, 'input/'t.bin' );
    out_file = 't.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( t, nx*ny*nz, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;

%% ------------------------------------------------------------------------
% Read out the noisy intital temperature to binary file.

if output;
    %out_file = fullfile( name.directory, name.problem, 'input/noisyt.bin' );
    out_file = 'noisyt.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( noisyt, nx*ny*nz, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;

%% ------------------------------------------------------------------------
% Read out the SURFACE FIELD of the noisy intital temperature to binary file.

if output;
    %out_file = fullfile( name.directory, name.problem, 'input/noisyt.bin' );
    out_file = 'noisytSurf.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( noisyt( :, :, 1 ), nx*ny, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;
%% ------------------------------------------------------------------------
% Read out the wind stress to binary file.

if output;
    out_file = fullfile( name.directory, name.problem, 'input/tau.bin' );
    out_file = 'tau.bin'
    % Open the file for writing only.
    [ fid message ] = fopen( out_file, 'w', format );
    % Reshape the bathymetry into a single column.
    output_data = reshape( tau, nx*ny, 1 );
    %Write it out.
    fwrite( fid, output_data, 'float32' );
    % Close the file.
    fclose( fid );
end;

%% ------------------------------------------------------------------------
% % Read out the surface heat flux to binary file.
% 
% if output;
%     out_file = fullfile( name.directory, name.problem, 'input/hflux.bin' );
%     % Open the file for writing only.
%     [ fid message ] = fopen( out_file, 'w', format );
%     % Reshape the bathymetry into a single column.
%     output_data = reshape( hflux, nx*ny, 1 );
%     %Write it out.
%     fwrite( fid, output_data, 'float32' );
%     % Close the file.
%     fclose( fid );
% end;

%% ------------------------------------------------------------------------
