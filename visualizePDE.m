function visualizePDE(ysol1, ysol2, t, nt)
%% this defines the spatial points for plotting purposes
% and is used in the pde solve
L = 1;
x = linspace(0,50,101);
t = linspace(0,t,nt);

% time step size used for plot titles
dx = t(2)-t(1);


%% this is were we solve for the entire solution
% you can replace the right hand sides with your solve
% sol1 = pdepe(m,@heatpde,@heatic,@heatbc1,x,t);
% sol2 = pdepe(m,@heatpde,@heatic,@heatbc2,x,t);
sol1 = ysol1;
sol2 = ysol2;

%% initializes the figure and size
figure(1)
set(gcf, 'Position',  [100, 100, 2500, 1000]);

%% the first two subplots are the entire time series of the solution
subplot(3,2,1)
colormap jet
imagesc(x,t,sol1)
% caxis([0 1])
colorbar
xlabel('Distance x','interpreter','latex')
ylabel('Time t','interpreter','latex')
title('Heat Equation for $0 \le x \le 1$ and $0 \le t \le 5$','interpreter','latex')

subplot(3,2,2)
colormap jet
imagesc(x,t,sol2)
% caxis([0 1])
colorbar
xlabel('Distance x','interpreter','latex')
ylabel('Time t','interpreter','latex')
title('Heat Equation for $0 \le x \le 1$ and $0 \le t \le 5$','interpreter','latex')

%% opens a video file for writing to
movie_name = 'solution';
vidfile = VideoWriter(sprintf('%s.mp4',movie_name),'MPEG-4');
open(vidfile);

%% iterate through the solution every '80' time steps in this case
% you can change the 80 to a higher number so that the video is shorter
% or you can decrease the number which will make the video longer
for i =1:5:length(t)
    if mod(i,1) == 0
        
        % this plots a 2d plot of position vs solution u
        s3=subplot(3,2,3);
        plot(x,sol1(i,:))
        xlabel('x')
        ylabel('u')
        %ylim([0 1])
        
        % this is a visualization of the movement of the solution u
        % along the 1 dimensional rod
        s5=subplot(3,2,5);
        scatter(x,x.*0,60,sol1(i,:),'filled','o')
        title(sprintf('t = %0.2f',i*dx));
        %caxis([0 1])
        colormap jet
        axis off
        
        % same as s3 subplot except for the second solution sol2
        s4 = subplot(3,2,4);
        plot(x,sol2(i,:))
        xlabel('x')
        ylabel('u')
        %ylim([0 1])
        
        s6 = subplot(3,2,6);
        scatter(x,x.*0,60,sol2(i,:),'filled','o')
        title(sprintf('t = %0.2f',i*dx));
        %caxis([0 1])
        colormap jet
        axis off
        
        drawnow
        
        % this captures the image frame to write to video files
        thisFrame = getframe(gcf);
        
        % write the image frame to the video file
        
        writeVideo(vidfile, thisFrame);
        
        % clear all axis
        cla(s3); cla(s5); cla(s4); cla(s6);
    end
end

% don't forget to close the video file!
close(vidfile);
end

%% these the pde functions you will need to replace these with your 
% calcium dynamic equations
function [c,f,s] = heatpde(x,t,u,dudx)
c = 100;
f = dudx;
s = 0;
end

function u0 = heatic(x)
u0 = 0.5;
end

function [pl,ql,pr,qr] = heatbc1(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur - 1;
qr = 0;
end

function [pl,ql,pr,qr] = heatbc2(xl,ul,xr,ur,t)
pl = ul;
ql = 0.75;
pr = ur - 1;
qr = 0.0;
end