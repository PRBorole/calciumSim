%function = calciumSim(filename,t,dt)

%{ 
 Relevant values from Breit, M., Kessler, M., Stepniewski, M. et al. 
 Spine-to-Dendrite Calcium Modeling Discloses Relevance for Precise Positioning of Ryanodine Receptor-Containing Spine Endoplasmic Reticulum. 
 Sci Rep 8, 15624 (2018). https://doi.org/10.1038/s41598-018-33343-9 
%}

clear all;
tic;
t = 30000; %simulation time
filename = "spatialresolvepde.swc";
stencil = stencilMakerCalcium(filename);
[~,~,~,~,R,~]=readSWC(filename); %Radii of Dendrite
[~,~,~,~,nx,~,meanEdge,~,~,~]=getGraphStructure(filename,false,false);
meanEdge = meanEdge*10^-5;
%nTrig = 10;
%ptLst = randi([1 nx],1,nTrig);
ptLst = [1];

R = R*10^-5;
r = R*15/40;

iv = getInitialValues();
nt = t/iv.dt; %Number of points in time discretization
iv.nx = nx;
iv.prob = iv.getInitProb(iv.Ceq);

% RyR initial probabilites
o1(1:nx,1) = iv.prob(1);
o2(1:nx,1) = iv.prob(2);
c1(1:nx,1) = iv.prob(3);
c2(1:nx,1) = 10^6-o1-o2-c1;

Ccy = iv.Ceq + (iv.Ceq*iv.kbPos*iv.b)/(iv.kbNeg + iv.Ceq*iv.kbPos); % Initial cytosolic Ca conc. 
Cey = iv.Ce + (iv.Ce*10^-5*iv.be)/(2000*10^-4 + iv.Ce*10^-5); % Initial ER Ca conc. 

Cc(1:nx,1) =  Ccy;%50*10^-18; %micromole/micrometer^3
b(1:nx,1) = iv.b;%40*10^-15; %micromole/micrometer^3, total CalB concentration, btot;
Ce(1:nx,1) = Cey;%250*10^-15; %micromole/micrometer^3
p(1:nx,1) = iv.p;%40*10^-18; %micromole/micrometer^3
be(1:nx,1) = iv.be;

RyRvalues = [o1 o2 c1 c2];
tSol = linspace(0,t,nt+1);
ySol = zeros(1,14*nx);

count = 1; 
for tstep = 1:iv.dt:t
    if mod(tstep,1000) == 0
        tstep
    end
    if tstep == 999
        [iv.vle,iv.vlp] = iv.getvlevlp(Cc(1),Ce(1),iv.Co,p(1),iv.prob');
        iv.Ce = Ce(1);
    end
    
    [Cc,b,Ce,p,be,RyRvalues,jP,jN,jR,jI,jS,jle,jlp,jERM,jSOC,jPM] = coupledEquations(tSol(count+1),tSol(count),[Cc b Ce p be],RyRvalues,stencil,R,r,meanEdge,ptLst,iv);
    ySol(count,:) = [Cc;b;Ce;p;be;jR;jI;jS;jle;jlp;jSOC;jPM;jP;jN];
    count = count + 1;
end

toc

tSol = tSol(1:nt);
list = ["Wave","Cc","b","Ce","be"];
plotFig(tSol,ySol,nx,iv.dt,list);
