% function [] = calciumSimTest()

%{ 
 Global Parameters

 Values from Breit, M., Kessler, M., Stepniewski, M. et al. 
 Spine-to-Dendrite Calcium Modeling Discloses Relevance for Precise Positioning of Ryanodine Receptor-Containing Spine Endoplasmic Reticulum. 
 Sci Rep 8, 15624 (2018). https://doi.org/10.1038/s41598-018-33343-9 
%}

global Dc Ccl Ccr;
Dc = 220; %micrometer^2/s, Diffusion constant
Ccl = 50; %nM, Left boundary
Ccr = 0; %nM, Right boundary
dt = 1*10^0; %time step
dx = 1; %x step
L = 50; %micrometer, lenght of neuron
t = 1000*10^-3; %second, simunlation time
nx = L/dx-1; %Number of points in space discretization
nt = t/dt; %Number of points in time discretization
alpha = Dc*dt/dx^2;

% Stencil with dirichlet's boundary condition
% stencil = zeros(nx+2);
% stencil(1,1) = 1;
% stencil(nx+2,nx+2) = 1;

% Stencil with boundary condition 
stencil = zeros(nx);
stencil(1,1) = 1+2*alpha;
stencil(1,2) = -alpha;
stencil(nx,nx) = -alpha;
stencil(nx,nx) = 1+2*alpha;

for i=2:nx-1
   for j=1:nx
       if i==j-1
           stencil(i,j) = -alpha;
       elseif i==j
           stencil(i,j) = 1+2*alpha;
       elseif i==j+1
           stencil(i,j) = -alpha;
       end
   end
end

% b, add matrix
b = zeros(nx,1);
add = zeros(nx,1);
add(1) = alpha*Ccl;
add(nx) = alpha*Ccr;


for i=1:nt
    b = inv(stencil)*(b+add);
end

plot(b)




