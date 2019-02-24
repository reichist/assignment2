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

% a) solve the simple case where V(x = 0) = V0 and V(x = L) = 0. Note that
% in this case the top/bottom BS are not fixed. You could use dV/dy = 0 for
% that BS or treat this as a 1D case. 

% Variables
W = 100e-9; 
L = 150e-9;
nx = (L*1e9)/2; 
ny = (W*1e9)/2;

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
            B(n) = 1;
        elseif z == nx
            G(n, n) = 1;
            B(n) = 0; 
        elseif p == 1
            G(n, n) = -3;
            G(n, nxm) = 1;
            G(n, nxp) = 1;
            G(n, nyp) = 1;
        elseif p == ny
            G(n, n) = -3;
            G(n, nxm) = 1;
            G(n, nxp) = 1;  
            G(n, nym) = 1;
        else
            G(n, n) = -4;
            G(n, nxm) = 1;
            G(n, nxp) = 1;
            G(n, nym) = 1;
            G(n, nyp) = 1;
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
title('Voltage Distribution - 1D');
xlabel('Length (u.a.)'); ylabel('Width (u.a.)'); zlabel('Voltage (V/V0)');
colorbar
view(135, 45)
