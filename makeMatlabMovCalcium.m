function makeMatlabMovCalcium(ySol, geometryfile,movieName)

% get coordinates of the geometry, the vertices 
[~,~,~,coords,r,~]=readSWC(geometryfile);

% Load time data from MatLab simulation
%filename=sprintf('%s/time.mat',dataFolder);

% get all the times for showing in plot
%load(filename,'t'); t = t';
tEnd = size(ySol,1);
nx = size(r,1);
t = linspace(0,tEnd,tEnd+1);

% this will make the geometry appear realistic by demonstrating the radii
% are not uniform, in particular the cell is thicker near the soma and
% thinner at the ends of dendrites
markerSize = (r./max(r)).*50;

% make a new figure windows
fig=figure('units','normalized','outerposition',[0 0 0.325 1.0]);
% this is for recording the movie
v = VideoWriter(sprintf('%s.mp4',movieName),'MPEG-4');
open(v)

% this part is for setting the text values on the color bar, no need to
% modify this except maybe the cmax and cmin if your action potentials have
% lower/higher peak values
yticklabel={};
cmax = 250; cmin = 0;
vals = [cmin:50:cmax];
for i=1:length(vals)
    yticklabel{i}=num2str(vals(i));
end
% if you change the 100 this will affect the length of the movie
for i=1:1:length(t)
        % read the voltage data from the .dat files
        u_sol = ySol(i,2*nx+1:3*nx);
        
        % make a scatter plot
        scatter3(coords(:,1),coords(:,2),coords(:,3),markerSize,'filled','CData',u_sol);
        
        %set labels
        xlabel(sprintf('{\\mu}m'))
        ylabel(sprintf('{\\mu}m'))
        set(gca,'Color', [0.5 0.5 0.5 0.5])
        caxis([cmin cmax])
        title(sprintf('MatLab, t = %0.2f [ms]',t(i+1)/10))
        colormap('jet')
        colorbar
        
        % set tick labels on colorbar
        c = colorbar;  
        c.Label.String="[microMolar]";
        c.YTick = [0:50:250];
        c.TickLabels = yticklabel; 
        view(2)
        
        % save the frame to video file
        thisframe=getframe(fig);
        writeVideo(v, thisframe);

        drawnow
        fprintf('frame = %i\n',i)
end
% don't forget to close the video file, if not it will be corrupted!
close(v)