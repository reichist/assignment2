clear all
clc

% Use the FD method to solve for the current flow in the rectangular region
% LxW shown in Fig 3 using del(sig_xy*delV) = 0. To model the boxes as
% highly resistive, we will use a low sig_xy inside the boxes. Note that
% sig_xy needs to remain finite. Start with sig = 1 outside the boxes and
% sig = 0.01 inside. In this exercise we'll experiment to see how the
% current flow is affected by the "bottle-neck". 

% (a) Calculate the current flow at the two contacts. Generate plots of
% sig(x, y), V(x, y), Ex, Ey, J(x, y). 

L = 150;
W = 100;

nx = L + 1;
ny = W + 1;
% nx = 2*L+1; %Investigate meshing
% ny = 2*W+1; %Investigate meshing

x = linspace(0, L, nx);
y = linspace(0, W, ny);

sigma = ones(nx, ny);
% sigma_in = 0.01;
sigma_in = 0.001; %investigate varying sigma in
% sigma_in = 0.1; %investigate varying sigma in
sigma_out = 1;

% Boxes Setup
box1coords = [((2*nx/6)) ((4*nx/6)) 0 (2*ny/6)];
box2coords = [((2*nx/6)) ((4*nx/6)) (4*ny/6) ny];
% box1coords = [((2*nx/6)) ((4*nx/6)) 0 (2.5*ny/6)]; %investigate bottleneck
% box2coords = [((2*nx/6)) ((4*nx/6)) (3.5*ny/6) ny]; %investigate bottleneck

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
            sigma(z, p) = sigma_out;
            G(n, n) = 1;
            B(n) = 1;
        elseif z == nx
            sigma(z, p) = sigma_out;
            G(n, n) = 1;
            B(n) = 1;
        elseif p == 1 || p == nx
            if z > box1coords(1) && z < box2coords(2)
                sigma(z, p) = sigma_in;
                G(n, n) = -3*sigma_in;
                G(n, nxm) = sigma_in;
                G(n, nxp) = sigma_in;
                G(n, nyp) = sigma_in;
            else
                sigma(z, p) = sigma_out;
                G(n, n) = -3*sigma_out;
                G(n, nxm) = sigma_out;
                G(n, nxp) = sigma_out;
                G(n, nyp) = sigma_out;
            end
        else
            G(n, n) = -4;
            if ((z > box1coords(1) && z < box1coords(2)) && (p > box1coords(3) && p < box1coords(4)))
                sigma(z, p) = sigma_in;
                G(n, nxm) = sigma_in;
                G(n, nxp) = sigma_in;
                G(n, nym) = sigma_in;
                G(n, nyp) = sigma_in;
            elseif ((z > box2coords(1) && z < box2coords(2)) && (p > box2coords(3) && p < box2coords(4))) 
                sigma(z, p) = sigma_in;
                G(n, nxm) = sigma_in;
                G(n, nxp) = sigma_in;
                G(n, nym) = sigma_in;
                G(n, nyp) = sigma_in;
            else
                sigma(z, p) = sigma_out;
                G(n, nxp) = sigma_out;
                G(n, nxm) = sigma_out;
                G(n, nyp) = sigma_out;
                G(n, nym) = sigma_out; 
            end
        end
    end
end

X = G\B;

voltage_distribution = zeros(nx, ny);
for z = 1:nx
    for p = 1:ny
        n = fn(z, p);
        voltage_distribution(z, p) = sigma(z, p) * X(n);
    end
end

figure(1)
clf
surf(sigma);
title('Sigma Matrix');
xlabel('nx'); ylabel('ny');
colorbar;


figure(2)
clf
surf(voltage_distribution);
title('Voltage Distribution (V/V0)');
xlabel('nx'); ylabel('ny');
view(135, 45);
colorbar;

[Ex, Ey] = gradient(voltage_distribution);

figure(3)
clf
subplot(1, 2, 1)
surf(Ex)
title('Ex');
xlabel('nx'); ylabel('ny');
colorbar;

% figure(5)
subplot(1, 2, 2)
surf(Ey)
title('Ey');
xlabel('nx'); ylabel('ny');
colorbar;

J = sigma.*gradient(voltage_distribution);

figure(4)
clf
surf(J)
title('Current Density');
xlabel('nx'); ylabel('ny');
colorbar;


% (b) Investigate mesh density



% (c) Investigate narrowing of "bottle-neck". 



% (d) Investigate varying of sigma inside the box. 