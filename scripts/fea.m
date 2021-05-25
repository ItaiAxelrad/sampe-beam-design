clc
clear all
close all

L = 23; %total length of highway overpass in inches
a = 9.5; if a >= L, error('a must be less than L'), end; %length from edge support to closest support in inches
b = L - a; if L - b ~= a, error('beam is not symmetrical'), end; %length from edge support to third support in feet
w = 750; %distributed loading on the overpass due to concrete and traffic in pounds per inch
E_Bamboo = 29000; %Modulus of elasticity for bamboo in psi
Y_Bamboo = 58; %Yielding strength for bamboo in psi
FS = 1.2; %Design factor of safety

A = [1 1 1 1; 0 a b L; 0 ((a^2) * (b^2)) / (3 * L) -(2 * (a^4) - (a^2) * (L^2)) / (6 * L) 0; ...
        0 -(((a^3) * L) - (3 * b * L * a^2) + (3 * a * L * b^2) + (2 * b^4) - (L * b^3) - ((b^2) * (L^2))) / (6 * L) ((a^2) * (b^2)) / (3 * L) 0];
B = [w * L; ((w * L^2) / 2); (w / 24) * ((a^4) - (2 * L * a^3) + (a * L^3)); (w / 24) * ((b^4) - (2 * L * b^3) + (b * L^3))];
R = LUfact(A, B); %R = reactions.

% Calculate and plot Shear and Moment
[min_shear, min_moment, max_shear, max_moment, V1func, V2func, ...
        V3func, M1func, M2func, M3func] = shear_moment(w, L, a, R(2), R(1));

% Find locations where moment is zero
root = zeros([1 6]); % Preallocate an array for roots
root(1) = fzero(M1func, [0 a / 2]);
root(2) = fzero(M1func, [a / 2 a]);
root(3) = fzero(M2func, [a L / 2]);
root(4) = fzero(M2func, [L / 2 L - a]);
root(5) = fzero(M3func, [L - a L - a / 2]);
root(6) = fzero(M3func, [L - a / 2 L]);
fprintf('\nMoment is zero at points (feet):\n')
fprintf('%10.3f  %10.3f  %10.3f %10.3f %10.3f %10.3f\n', root)

% Perform analyses on structural steel
fprintf('\n Structural Steel \n')
[name, A, d, bf, tf, tw, Ixx] = wbeamselect(max(abs([min_moment max_moment])), ...
    Y_Steel, FS);

% Calculate and plot slope and displacement of the beam
[min_slope, min_displacement, max_slope, max_displacement] = ang_disp(w, R(2), L, a, E_steel, Ixx, 'Structural Steel');

% Calculate stress at arbitrary state (100ft,9in) and angle theta of 45 degrees
theta_d = 45;
[sigma_xxp, sigma_yyp, tau_xyp, min_princ, max_princ, thetap1, thetap2] ...
    = stressstate(100, 9, V1func, M1func, [0 a], V2func, M2func, [a b], V3func, M3func, [b L], ...
    theta_d, d, bf, tf, tw, Ixx);
fprintf('\n sigma xx (ksi)    sigma yy(ksi)   tau xy(ksi)   angle(degrees)\n')
fprintf('%10.3f    %10.3f      %10.3f     %10.3f\n', sigma_xxp, sigma_yyp, tau_xyp, theta_d)
fprintf('\n Principle Stresses (ksi)         Principle Planes (degrees) \n')
fprintf('%10.3f     %10.3f      %10.3f     %10.3f\n', min_princ, max_princ, thetap1, thetap2)

% Plot stresses with respect to angle and Mohr's Circle at (100ft,9in)
transplot(100, 9, V1func, M1func, [0 a], V2func, M2func, [a b], V3func, M3func, [b L], ...
    d, bf, tf, tw, Ixx, 'Structural Steel');

% Make false-color plot of Factor of Safety
FS = fsplot(Y_Steel, V1func, M1func, [0 a], V2func, M2func, [a b], V3func, M3func, [b L], ...
    d, bf, tf, tw, Ixx, 'Structural Steel');
fprintf('\n Overall Factor of Safety: %f\n', FS)
