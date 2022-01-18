% function [] = calciumSimTest()

%{ 
 Global Parameters

 Relevant values from Breit, M., Kessler, M., Stepniewski, M. et al. 
 Spine-to-Dendrite Calcium Modeling Discloses Relevance for Precise Positioning of Ryanodine Receptor-Containing Spine Endoplasmic Reticulum. 
 Sci Rep 8, 15624 (2018). https://doi.org/10.1038/s41598-018-33343-9 
%}

clear all;
tic;
dt = 1; %time step 0.0625
t = 10000; %msecond, simulation time
nt = t/dt; %Number of points in time discretization
% stencil = stencilMakerCalcium("sample_geometry/Cells_for_Piyush/Cell1/228-16MG.CNG_segLength=8_1d_ref_2.swc");
filename = "spatialresolvepde.swc";
% filename = "sample_geometry/Cells_for_Piyush/Cell7/194-4-17nj.CNG_segLength=12_1d_ref_2.swc";
% filename = "yShaped1.swc";
stencil = stencilMakerCalcium(filename);
% stencil = stencilMakerCalcium("2input2output.swc");
[~,~,~,~,R,~]=readSWC(filename); %Radii of Dendrite
[~,~,~,~,nx,~,meanEdge,~,~,~]=getGraphStructure(filename,false,false);
meanEdge = meanEdge*10^-5;
nTrig = 10;
%ptLst = randi([1 nx],1,nTrig);
ptLst = [1];

% idx = randi([1 5], 1,1);
% ptLst = ptLst(idx);

R = R*10^-5;
r = R*15/40;
% R = 4*R*10^-5; % radii of Dendrite
% r = R*15/40; % radii of ER
% r(1:151) = 0.02*10^-5;
% r(202:351) = 0.02*10^-5;

%RyR intital probabilites
o1(1:nx,1) = 324;
o2(1:nx,1) = 0;
c1(1:nx,1) = 994014;
c2(1:nx,1) = 10^6-o1-o2-c1;

Ceq = 0.05; %mM
kbNeg = 19*10^-4;%19; %sec^-1; CalB used 
kbPos = 27*10^-4;%27*10^15; %micrometer^3/micromol s; CalB free
Ccy = Ceq + (Ceq*kbPos*40)/(kbNeg + Ceq*kbPos);

Cc(1:nx,1) =  Ccy;%50*10^-18; %micromole/micrometer^3
b(1:nx,1) = 40;%40*10^-15; %micromole/micrometer^3, total CalB concentration, btot;
Ce(1:nx,1) = 250;%250*10^-15; %micromole/micrometer^3
p(1:nx,1) = 40*10^-3;%40*10^-18; %micromole/micrometer^3

RyRvalues = [o1 o2 c1 c2];
tSol = linspace(0,t,nt+1);
ySol = zeros(1,7*nx);

count = 1;
for tstep = 1:dt:t
    if mod(tstep,1000) == 0
        tstep
    end
    if tstep == 1000
        %RyR intital probabilites
        o1(1:nx,1) = 324;
        o2(1:nx,1) = 0;
        c1(1:nx,1) = 994014;
        c2(1:nx,1) = 10^6-o1-o2-c1;
        RyRvalues = [o1 o2 c1 c2];
    end
    
    [Cc,b,Ce,p,RyRvalues,jP,jN,jR,jI,jS,poI] = coupledEquations(tSol(count+1),tSol(count),[Cc b Ce p],nx,RyRvalues,stencil,R,r,meanEdge,ptLst);
    ySol(count,:) = [Cc;b;Ce;p;jR;jI;poI];
    count = count + 1;
end

toc

tSol = tSol(1:nt);
% figure;
% plot(tSol,ySol(:,1:nx))
% title("ODE45-Cc");
% % legend();
% 
% figure;
% plot(tSol,ySol(:,nx+1:2*nx))
% title("ODE45-b");
% % legend();
% 
% figure;
% plot(tSol,ySol(:,2*nx+1:3*nx))
% title("ODE45-Ce");
% % legend();
% 
% figure;
% plot(tSol,ySol(:,3*nx+1:4*nx))
% title("ODE45-p");
% legend();


figure;
plot(tSol,ySol(:,1))
title("ODE45-Cc");
% legend();

figure;
plot(tSol,ySol(:,nx+1))
title("ODE45-b");
% legend();

figure;
plot(tSol,ySol(:,2*nx+1))
title("ODE45-Ce");
% legend();

figure;
plot(tSol,ySol(:,3*nx+1))
title("ODE45-p");
legend();

% figure; plot(linspace(0,64,nx),ySol(1500/dt:10:t/dt,1:nx))
% figure;
% plot(linspace(0,100,201),ySol(1500/dt:100:25000/dt,1:201))
% title("Cytosol Ca^{2+} vs Length 9pt");
% legend(string([1500/dt:10:3000/dt]));

% figure; plot(linspace(0,100,201),ySol(1500/dt:100:25000/dt,[351:-1:202,151:201]))

% figure;
% diff = pt(nx-64)-pt(nx-65);
% plot(linspace(0,diff*64,64),ySol(2900:50:6000,nx-63:nx))
% title("Cytosol Ca^{2+} vs Length 193pt");
% legend(string([2900:50:6000]));

% figure;
% plot(tSol,ySol(:,4*nx+1:5*nx))
% title("jI");
% legend();
% 
% figure;
% plot(tSol,ySol(:,5*nx+1:6*nx))
% title("ODE45-p");
% legend();
% 
% figure;
% plot(tSol,ySol(:,6*nx+1:7*nx))
% title("ODE45-jI");
% legend();
% 
% figure;
% plot(tSol,ySol(:,7*nx+1:8*nx))
% title("ODE45-jS");
% legend();






