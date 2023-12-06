%========================================================================
%
%   Model of free motion on a rotating frame
%                                          _
%   du              dv                      |
%   -- - fv = 0,    -- + fu = 0,    t > 0   |_
%   dt              dt                      |
%   u(0) = u_0,     v(0) = v_0,     t = 0  _|
%
%   + Traking                    _
%   dx         dy                 |
%   -- = u,    -- = v,    t > 0   |_
%   dt         dt                 |
%   x(0) = 0,  y(0) = 0,  t = 0  _|
%
% Copyright © 2023 Llorente Lázaro, Víctor Javier                         
% Last Update: Dec 6, 2023                                               
% Website: https://sites.google.com/view/vjllorente                       
% Contact: victor.javier.llorente@gmail.com  
%
%==========================================================================

clear all
close all
clc

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  INPUTS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial velocity field [m/s]
u0 = 1;
v0 = 1;
% Latitud [º ' '']
lat = 37 + 1 / 60 + 0 / 3600;
% Earth's period and rotation
T     = 23 * 3600 + 56 * 60 + 4.1;                                  
Omega = 2 * pi / T;       
% Coriolis factor [1/s]
f = 2 * Omega * sin( lat * pi / 180 );
% End time [s]
tend = 1000000;
% Loop time [s]
tloop = 2 * pi / f;
% Time step [s]
Dt = 100;
% Numerical scheme 
% exp    = explicit            (unstable)
% imp    = implicit            (overly stable)
% expimp = explicit + implicit (oscillatory, condition: fDt < 2)
% fis    = Fischer             (energy stable)
esquema = 'exp';
if strcmp( esquema, 'expimp' )
    Dt = 0.6 / f;
end
% Dimensionaless parameter
r = f * Dt;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                  MESH                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discrete time [s]
tn = [ 0 : Dt : tend ];
% Number of nodes
N = length( tn );

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               VELOCITY FIELD                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate
u = zeros( 1, N );
v = zeros( 1, N );
% Initial velocities
u( 1 ) = u0;
v( 1 ) = v0;
% Main loop
switch esquema
    case 'exp'
        for n = 2 : N
            u( n ) = u( n - 1 ) + r * v( n - 1 );
            v( n ) = v( n - 1 ) - r * u( n - 1 );
        end        
    case 'imp'
        for n = 2 : N
            u( n ) = ( u( n - 1 ) + r * v( n - 1 ) ) / ( 1 + r ^ 2 );
            v( n ) = ( v( n - 1 ) - r * u( n - 1 ) ) / ( 1 + r ^ 2 );
        end
    case 'expimp'
        for n = 2 : N
            u( n ) = u( n - 1 ) + r * v( n - 1 );
            v( n ) = v( n - 1 ) - r * u( n     );
        end   
    case 'fis'
        f1 = ( 1 - ( r / 2 ) ^ 2 ) / ( 1 + ( r / 2 ) ^ 2 );
        f2 = r                     / ( 1 + ( r / 2 ) ^ 2 );
        for n = 2 : N
            u( n ) = f1 * u( n - 1 ) + f2 * v( n - 1 );
            v( n ) = f1 * v( n - 1 ) - f2 * u( n - 1 );
        end
    otherwise
        warning('No correct numerical scheme')
end
% Kinetic energy
k = ( u .^ 2 + v .^ 2 ) ./ 2;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 TRAKING                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Allocate
x = zeros( 1, N );
y = zeros( 1, N );
% Start point
x( 1 ) = 0;
y( 1 ) = 0;
% Integration of lagrangian formulation by the Trapezoidal rule
for n = 2 : N 
    x( n ) = x( n - 1 ) + Dt * ( u( n ) + u( n - 1 ) ) / 2;
    y( n ) = y( n - 1 ) + Dt * ( v( n ) + v( n - 1 ) ) / 2;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              EXACT SOLUTIONS                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact velocity field
ue = u0 .* cos( f .* tn ) + v0 .* sin( f .* tn );
ve = v0 .* cos( f .* tn ) - u0 .* sin( f .* tn );
% Exact kinetic energy
ke = ( ue .^ 2 + ve .^ 2 ) ./ 2;
% Trayectory
t = [0:1000:tloop];
xe = ( u0 / f ) .* sin( f .* t ) + ( v0 / f ) .* ( 1 - cos( f .* t ) );
ye = ( v0 / f ) .* sin( f .* t ) - ( u0 / f ) .* ( 1 - cos( f .* t ) );
% Center and radius of curvature: ( xe - xo ) ^ 2 + ( ye - yo ) ^ 2 = R^2
xo =  v0 / f; 
yo = -u0 / f;
R = sqrt( u0^2 + v0^2 ) / f;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   PLOTS                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1: trayectory
figure( 1 )
plot( x , y , '-b' ); 
hold on
plot( xe, ye, '--r' );
legend( 'Numerical', 'Exact' )
xlabel( '$x$ [m]', 'interpreter', 'latex' )
ylabel( '$y$ [m]', 'interpreter', 'latex' )
grid on
set( gcf, 'color', 'w' )
% Figure 2: Kinetic energy vs time
figure( 2 )
plot( tn, k, '-b', tn, ke, '-r' );
grid on
legend( 'Numerical', 'Exact' )
xlabel( '$t$ [s]', 'interpreter', 'latex' )
ylabel( '$(u^2+v^2)/2$ [m$^2$/s$^2$]', 'interpreter', 'latex' )
grid on
set( gcf, 'color', 'w' )
% Figure 3: Velocity field vs time
figure( 3 )
subplot(2,1,1);
plot( tn, u, '-b', tn, ue, '-r' )
grid on
legend( 'Numerical', 'Exact' )
xlabel( '$t$ [s]', 'interpreter', 'latex' )
ylabel( '$u$ [m/s]', 'interpreter', 'latex' )
grid on
set( gcf, 'color', 'w' )
subplot(2,1,2); 
plot( tn, v, '-b', tn, ve, '-r' )
grid on
legend( 'Numerical', 'Exact' )
xlabel( '$t$ [s]', 'interpreter', 'latex' )
ylabel( '$v$ [m/s]', 'interpreter', 'latex' )
grid on
set( gcf, 'color', 'w' )
