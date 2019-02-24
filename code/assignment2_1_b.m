clear all
clc

% Laplace's equation by FD can be used to solve electrostatic potential
% problems del2V = 0, or current flow problems in inhomogeneous solids
% del(sig_xydelV) = 0. Either case can be modeled as an orthogonal resistor
% network with resistors of value 1 or sig^-1. The conception of the
% problem is an orthogonal 2D mesh resistors is very useful for setting up
% the boundary conditions (BC) and the interface conditions where 2 regions
% have different sigmas

% 1. Use the Finite Difference Method to solve for the electrostatic
% potential in the rectangular region LxW shown in Figure 1 using del2V = 0
% in this excercise we'll experiment to see how accurate our model is. 

% a) solve the case where V(x = 0) = V0 and V(x = L) = V0 and V(y = 0) = 0 
% and V(y = W) = 0. Note that in this case the top/bottom BS are not fixed.
% You could use dV/dy = 0 for that BS or treat this as a 1D case. 

% Variables
W = 100; 
L = 150;

nx = L/2;
ny = W/2;

V0 = 1;

lBC = V0;
rBC = V0;
uBC = 0; 
dBC = 0;

G = sparse(nx*ny, nx*ny);
B = zeros(nx*ny, 1);

fn = @(i, j) j + (i-1)*ny;

for z = 1:nx
    for p = 1:ny
        n = fn(z, p);
        nxm = fn(z-1, p);
        nxp = fn(z+1, p);
        nym = fn(z, p-1);
        nyp = fn(z, p+1);
        if z == 1
            G(n, n) = 1;
            B(n) = lBC;
        elseif z == nx
            G(n, n) = 1;
            B(n) = rBC;
        elseif p == 1
            G(n, n) = 1;
            B(n) = dBC;
        elseif p == ny
            G(n, n) = 1;
            B(n) = uBC;
        else 
            G(n, n) = -4;
            G(n, nxp) = 1;
            G(n, nxm) = 1;
            G(n, nyp) = 1;
            G(n, nym) = 1; 
        end
    end
end

X = G\B;
voltage_distribution = zeros(nx, ny);
for z = 1:nx
    for p = 1:ny
        n = fn(z, p);
        voltage_distribution(z, p) = X(n);
    end
end

figure(1)
clf
surf(voltage_distribution)
title('Voltage Distribution');
xlabel('Length (a.u.)'); ylabel('Width (a.u.)'); zlabel('Voltage (V/V0)');
colorbar;

% Compare the solution of a bunch of mesh sizes to the analytical series
% solution. Drwa some conclusions about meshing. Note that both solutions
% have some explicit error. It is useful to plot a "movie" of the
% analytical series solution as it converges to the "real" solution. How do
% you know/decide when to stop the analytical series? Griffiths "Intro to
% Electrodynamcis 3e" contains a derivation of the analytical solution in
% Example 3.4/ Comment on the advantages and disadvantages of the numerical
% versus analytical solution. 

% http://ed.quantum-bg.org/Introduction%20To%20Electrodynamics%20-Griffiths.pdf
% slide 141

Vconst = 4*V0/pi;
% n = 1;
% n = [1 3];
n1 = 1:2:10;
n2 = 1:2:302;
X1 = zeros(nx, ny);
X2 = zeros(nx, ny);
a = W;
b = L/2;
x_axis = linspace(-b, b, nx);
y_axis = linspace(0, a, ny);

Voltage1 = @(x, y, n) cosh((n*pi*x/a));
Voltage2 = @(x, y, n) cosh((n*pi*b/a));
Voltage3 = @(x, y, n) sin((n*pi*y/a));

v1111 = Voltage1(1, 1, 1);
v2111 = Voltage2(1, 1, 1);
v3111 = Voltage3(1, 1, 1);

Voltage = @(x, y, n) ((Voltage1(x, y, n)/Voltage2(x, y, n)*Voltage3(z, y, n)))/n;
v111 = Voltage(1, 1, 1);

for z = 1:length(x_axis)
    for p = 1:length(y_axis)
        sum1 = 0;
        sum2 = 0;
        for w = 1:length(n1)
            sum1 = sum1 + Voltage(x_axis(z), y_axis(p), n1(w));
        end
        for w = 1:length(n2)
            sum2 = sum2 + Voltage(x_axis(z), y_axis(p), n2(w));
        end
        X1(z, p) = sum1*Vconst;
        X2(z, p) = sum2*Vconst;
    end
    
    figure(2)
    subplot(1, 2, 1)
    surf(X1)
    title(strcat('(n=', num2str(length(n1)),' sums/node)'));
    xlabel('x (a.u.)'); ylabel('y (a.u)'); zlabel('Voltage (V/V0)');
    % view(2);
    
    figure(2)
    subplot(1, 2, 2)
    surf(X2)
    title(strcat('(n=', num2str(length(n2)),' sums/node)'));
    xlabel('x (a.u.)'); ylabel('y (a.u)'); zlabel('Voltage (V/V0)');
    % view(2);
    
    pause(0.01)    
end

figure(4)
clf
surf(X1)
title(strcat('Analytical Series Solution (n=', num2str(length(n1)),' sums/node)'));
xlabel('x (a.u.)'); ylabel('y (a.u)'); zlabel('Voltage (V/V0)');
colorbar;

figure(3)
clf
surf(X2)
title(strcat('Analytical Series Solution (n=', num2str(length(n2)),' sums/node)'));
xlabel('x (a.u.)'); ylabel('y (a.u)'); zlabel('Voltage (V/V0)');
colorbar;
    
figure(5)
clf
subplot(1, 2, 1)
surf(voltage_distribution);
title('GX=B Method');
xlabel('Length (a.u.)'); ylabel('Width (u.a.)'); zlabel('Voltage (V/V0)');
colorbar;

subplot(1, 2, 2)
surf(X2)
title(strcat('Analytical Series Solution (n=', num2str(length(n2)),' sums/node)'));
xlabel('x (a.u.)'); ylabel('y (a.u)'); zlabel('Voltage (V/V0)');
colorbar;
