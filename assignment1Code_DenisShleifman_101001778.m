%% Main Code

function [] = assignment1Code_DenisShleifman_101001778()
clear
clc
close all

const.me = 9.10938356e-31; %electron mass(kg)
const.re = 2.8179403227e-15; %electron radius (m)
const.kb = 1.38066e-23; %Boltzmans Constant (J/K)
const.T = 300; %Board Temperature (K)
const.tmn = 0.2e-12; % Mean time between collisions

sim.particles = 1000; %number of particles to initialize in the simulation
sim.delay = 0.01; %time in seconds between executing steps
sim.steps = 1000; %number of iterations to perform
sim.movie = true; %set to true if you would like it to loop through each iteration
sim.question = 3; %set which question on the assignment to simulate
sim.scatter = true; %set whether boundaries cause
sim.heatmapBins = [10,10]; %number of bins for 2D heatmap [x,y]

sim.dist.mult = 1e9;
sim.dist.unit = 'nm';
sim.time.mult = 1e12;
sim.time.unit = 'ps';
sim.speed.mult = 1e-3;
sim.speed.unit = 'km/s';

part.colour = blanks(sim.particles);
part.colour(1) = 'b';
part.colour(2) = 'g';
part.colour(3) = 'r';
part.colour(4) = 'c';
part.colour(5) = 'm';
part.colour(6) = 'y';
part.colour(7) = 'k';

board.xMin = 0; %Minimum X coordinate for the board (m)
board.yMin = 0; %Minimum Y coordinate for the board (m)
board.xMax = 200e-9; %Maximum X coordinate for the board (m)
board.yMax = 100e-9; %Maximum Y coordinate for the board (m)

performSimuation(const,part,board,sim);
end

function [] = performSimuation(const,part,board,sim)
[const,part,board,sim,stat,plt] = initialize(const,part,board,sim);

for n = 1:sim.steps
    [part,stat,plt,sim] = performIteration(const,part,board,sim,stat,plt,n);
end
updatePlot(plt,part,stat,sim,sim.steps,true);
outputCalcResults(const,sim,stat);
plotHeatMaps(board,sim,stat);
end
%% Initialization of all variables

function [const,part,board,sim,stat,plt] = initialize(const,part,board,sim)
board = initializeBoard(board,sim.question);
const = initializeConst(const);
sim = initializeSim(sim,const,board.minSideLength);
part = initializePart(part,board,sim,const);
stat = initializeStats(board,sim);
plt = setupPlot(board,sim,const);
end

function [board] = initializeBoard(board,question)
board.Lx = board.xMax - board.xMin;
board.Ly = board.yMax - board.yMin;
board.A = board.Lx*board.Ly;
board.minSideLength = min([board.Lx,board.Ly]);
board = initializeBoundaries(board,question);
end

function [board] = initializeBoundaries(board,question)
startX = board.xMin - board.Lx/5;
endX = board.xMax + board.Lx/5;

if question ~= 3
    board.boundaries{1} = [startX,board.yMin,endX,board.yMin];
    board.boundaries{2} = [startX,board.yMax,endX,board.yMax];
    return
end
xBoxLeft = board.xMin + board.Lx*2/5;
xBoxRight = board.xMin + board.Lx*3/5;
yBoxTop = board.yMin + board.Ly*4/9;
yBoxBot = board.yMin + board.Ly*5/9;

board.boxes{1} = [xBoxLeft,board.yMin,xBoxRight,yBoxTop];
board.boxes{2} = [xBoxLeft,yBoxBot,xBoxRight,board.yMax];
board.boundaries{1} = [startX,board.yMin,xBoxLeft,board.yMin];
board.boundaries{2} = [xBoxLeft,board.yMin,xBoxLeft,yBoxTop];
board.boundaries{3} = [xBoxLeft,yBoxTop,xBoxRight,yBoxTop];
board.boundaries{4} = [xBoxRight,yBoxTop,xBoxRight,board.yMin];
board.boundaries{5} = [xBoxRight,board.yMin,endX,board.yMin];
board.boundaries{6} = [startX,board.yMax,xBoxLeft,board.yMax];
board.boundaries{7} = [xBoxLeft,board.yMax,xBoxLeft,yBoxBot];
board.boundaries{8} = [xBoxLeft,yBoxBot,xBoxRight,yBoxBot];
board.boundaries{9} = [xBoxRight,yBoxBot,xBoxRight,board.yMax];
board.boundaries{10} = [xBoxRight,board.yMax,endX,board.yMax];
board.corners{1} = [xBoxLeft,board.yMin];
board.corners{2} = [xBoxRight,board.yMin];
board.corners{3} = [xBoxLeft,board.yMax];
board.corners{4} = [xBoxRight,board.yMax];

for n = 1:length(board.boxes)
    board.A = board.A - abs((board.boxes{n}(3) - board.boxes{n}(1))*(board.boxes{n}(4) - board.boxes{n}(2)));
end
end

function [const] = initializeConst(const)
const.meff = const.me*0.26;
const.de = const.re*2;
const.vth = sqrt((2*const.T*const.kb)/const.meff);
const.vav = const.vth*sqrt(pi/4);
const.MFP = const.vav *const.tmn;
end

function [sim] = initializeSim(sim,const,minSideLength)
sim.dt = minSideLength/(100*const.vth);
sim.outputPercent = 0;
sim.Pscat = 1 - exp(-sim.dt/const.tmn);
end

function [part] = initializePart(part,board,sim,const)
part = generateInitialPosition(part,board,sim);
part = assignVelocities(part,sim,const);
part.colTime = zeros(1,sim.particles);
part.colDist = zeros(1,sim.particles);
part.v = zeros(1,sim.particles);
part.v2 = zeros(1,sim.particles);
part.temp = zeros(1,sim.particles);
part.traceIndex = part.colour ~= ' ';
end

function [part] = assignVelocities(part,sim,const,I)
if nargin == 3
    I = true(1,sim.particles);
end

if sim.question == 1
    angle = generatePoints(sum(I),0,2*pi);
    part.vx(I) = cos(angle)*const.vth;
    part.vy(I) = sin(angle)*const.vth;
else
    part.vx(I) = nrmrnd(0,1,1,sum(I))*const.vth/sqrt(2);
    part.vy(I) = nrmrnd(0,1,1,sum(I))*const.vth/sqrt(2);
end
end

function [out] = generatePoints(N,Min,Max)
out = rand(1,N)*(Max - Min) + Min;
end

function nrmmatrix = nrmrnd(mu, sigma, sz1, sz2)
nrmmatrix = mu+sigma*randn(sz1,sz2);
end

function [part] = generateInitialPosition(part,board,sim)
I = true(1,sim.particles);
while sum(I) ~= 0
    part.x(I) = generatePoints(sum(I),board.xMin,board.xMax);
    part.y(I) = generatePoints(sum(I),board.yMin,board.yMax);
    I = false(1,sim.particles);
    if sim.question == 3
        for n = 1:length(board.boxes)
            I = I | isInsideBox(board.boxes{n},part.x,part.y);
        end
    end
end
end

function [I] = isInsideBox(box,x,y)
I = (x > box(1) & x < box(3) & y > box(2) & y < box(4));
end

function [stat] = initializeStats(board,sim)
stat.time = zeros(1,sim.steps);
stat.temp = zeros(1,sim.steps);
stat.avtemp = zeros(1,sim.steps);
stat.MFP = zeros(1,sim.steps);
stat.tmn = zeros(1,sim.steps);
stat.avV = zeros(1,sim.steps);
stat.avVrms = zeros(1,sim.steps);
stat.heatMap = generateHeatMapXY(board,sim);
stat.avHeatMap = zeros(sim.heatmapBins(1)+1,sim.heatmapBins(2)+1);
stat.avTempMap = zeros(sim.heatmapBins(1)+1,sim.heatmapBins(2)+1);
end

function [heatMap] = generateHeatMapXY(board,sim)
heatMap.x = linspace(board.xMin,board.xMax,sim.heatmapBins(1)+1);
heatMap.y = linspace(board.yMin,board.yMax,sim.heatmapBins(2)+1);
[heatMap.X,heatMap.Y] = meshgrid(heatMap.x*sim.dist.mult,heatMap.y*sim.dist.mult);
end

%% Setup Plots

function [plt] = setupPlot(board,sim,const)
plt.fig1 = setupParticlePlot(board,sim);
if sim.question ~= 3
    plt.fig2 = setupTimePlot(sim);
end
if sim.question == 2
    plt.fig3 = setupVelocityDist(sim,const);
end
end

function [fig] = setupParticlePlot(board,sim)
fig.h = figure;
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
axis(fig.ax1,[board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
hold(fig.ax1, 'on');
if sim.movie
    fig.plotPos = plot(fig.ax1,1,1,'.b');
    fig.plotVel = quiver(fig.ax1,[],[],[],[],2,'b');
end
drawBoundaries(fig.ax1,sim,board.boundaries);
xlabel(fig.ax1,sprintf('x (%s)',sim.dist.unit));
ylabel(fig.ax1,sprintf('y (%s)',sim.dist.unit));
title(fig.ax1,sprintf('Particle Trajectories (7/%d)',sim.particles));
grid(fig.ax1,'on');
end

function [] = drawBoundaries(ax,sim,boundaries)
if sim.question == 3
    for n = 1:length(boundaries)
        currentBoundary = boundaries{n};
        line(ax,[currentBoundary(1),currentBoundary(3)]*sim.dist.mult,[currentBoundary(2),currentBoundary(4)]*sim.dist.mult,'Color','k');
    end
end
end

function [fig] = setupTimePlot(sim)
fig.h = figure();
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
hold(fig.ax1, 'on')
fig.plotTemp = plot(fig.ax1,1,1,'-b');
if sim.question ~= 1
    fig.plotAvTemp = plot(fig.ax1,1,1,'-r');
    legend(fig.ax1,{'Current','Average'});
end
xlabel(fig.ax1,sprintf('Time (%s)',sim.time.unit));
ylabel(fig.ax1,'Temperature (K)');
grid(fig.ax1,'on');
end

function [fig] = setupVelocityDist(sim,const)
fig.h = figure();
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
fig.axis = axis(fig.ax1);
hold(fig.ax1, 'on')
fig.vHist = histogram(fig.ax1,[],'FaceColor','g');
xline(fig.ax1,const.vth*sim.speed.mult,'Label','vth','Color','r','LineWidth',2);
fig.Vav = xline(fig.ax1,0,'Label','Average Velocity','Color','b','LineWidth',2,'LabelVerticalAlignment','bottom');
fig.Vrms = xline(fig.ax1,0,'Label','Vrms','Color','k','LineWidth',2,'LabelVerticalAlignment','middle');
xlabel(fig.ax1,sprintf('Velocity (%s)',sim.speed.unit));
ylabel(fig.ax1,'Counts');
grid(fig.ax1,'on');
end

%%Perform Iteration

function [part,stat,plt,sim] = performIteration(const,part,board,sim,stat,plt,n)
part = advancePosition(part,board,sim,const);
stat = updateStats(stat,part,sim.dt,n);
plt = updatePlot(plt,part,stat,sim,n);
if sim.movie
    pause(sim.delay);
else
    sim = outputPercentDone(sim,n);
end
end

function [sim] = outputPercentDone(sim,currStep)
while (currStep/sim.steps)*100 > sim.outputPercent
    sim.outputPercent = sim.outputPercent + 1;
    fprintf('Percent Complete: %d %%\n',sim.outputPercent);
end
end

function [part] = advancePosition(part,board,sim,const)
lastX = part.x;
lastY = part.y;

part = boundaryCondition(part,board,sim,const);
part = calcValues(part,const,sim.dt);
part = updateAdd(part,lastX,lastY,board.Lx);
part = scatter(part,sim,const);
end

function [part] = boundaryCondition(part,board,sim,const)
dt = ones(1,sim.particles) * sim.dt;
forbid = true(1,sim.particles);
while sum(forbid) ~= 0
    [x,y] = calcNewPos(part,dt);
    forbid = checkForbiddenRedions(sim,board,x,y);
    part.x(~forbid) = x(~forbid);
    part.y(~forbid) = y(~forbid);
    dt(~forbid) = 0;
    if sim.question == 3
        for n = 1:length(board.corners)
            [part,dt] = checkCorners(part,sim,board,const,dt,forbid,board.corners{n},x,y);
        end
    end
    for n = 1:length(board.boundaries)
        [part,dt] = checkBoundaries(part,sim,const,dt,forbid,board.boundaries{n},x,y);
    end
    forbid = checkForbiddenRedions(sim,board,x,y);
end
part.x = wrap(part.x,board.xMin,board.xMax);
end

function [x,y] = calcNewPos(part,dt)
x = part.x + part.vx.*dt;
y = part.y + part.vy.*dt;
end

function [forbid,x,y] = checkForbiddenRedions(sim,board,x,y)
forbid = y > board.yMax | y < board.yMin;
if sim.question == 3
    for n = 1:length(board.boxes)
        forbid = forbid | isInsideBox(board.boxes{n},x,y);
    end
end
end

function [part,dt] = checkCorners(part,sim,board,const,dt,forbid,currentCorner,x,y)
intersect = pointIntersection(currentCorner,[part.x;part.y;x;y],forbid);
if sum(intersect) == 0
    return;
end
dist2 = (currentCorner(1) - part.x(intersect)).^2 + (currentCorner(2) - part.y(intersect)).^2;
v2 = part.vx(intersect).^2 + part.vy(intersect).^2;
dt(intersect) = dt(intersect) - sqrt(dist2./v2);
part.x(intersect) = currentCorner(1);
part.y(intersect) = currentCorner(2);
allowedReflection = ~checkForbiddenRedions(sim,board,currentCorner(1) + board.Lx*[1,-1,-1,1]/100,currentCorner(2) + board.Ly*[1,1,-1,-1]/100);
if sim.question == 3 && sim.scatter
    part = assignScatterBoundaryVelCorner(part,sim,const,intersect);
else
    if allowedReflection(1) || allowedReflection(3) %corner is region 1 or region 3
        part.vx(intersect) = -part.vy(intersect);
        part.vy(intersect) = -part.vx(intersect);
    else
        part.vx(intersect) = part.vy(intersect);
        part.vy(intersect) = part.vx(intersect);
    end
end
end

function [part] = assignScatterBoundaryVelCorner(part,sim,const,intersect)
expectedSignsX = zeros(1,sim.particles);
expectedSignsY = zeros(1,sim.particles);
expectedSignsX(intersect) = -sign(part.vx(intersect));
expectedSignsY(intersect) = -sign(part.vy(intersect));
IAssign = intersect;
while sum(IAssign) ~= 0
    part = assignVelocities(part,sim,const,IAssign);
    IAssign(intersect) = sign(part.vx(intersect)) ~= expectedSignsX(intersect) & sign(part.vy(intersect)) ~= expectedSignsY(intersect);
end
end

function [intersect] = pointIntersection(point,lines,include)
if nargin == 2
    include = true(1,size(lines,2));
end
intersect = false(1,length(include));

dx = lines(3,include) - lines(1,include);
dy = lines(4,include) - lines(2,include);
dy2 = point(2) - lines(2,include);
dx2 = point(1) - lines(1,include);

intersect(include) = dy.*dx2 == dy2.*dx & point(2) >= lines(1,include) & point(2) <= lines(3,include);
end

function [part,dt] = movePointToCrossing(part,dt,intersect,ix,iy)
dist2 = (ix - part.x(intersect)).^2 + (iy - part.y(intersect)).^2;
v2 = part.vx(intersect).^2 + part.vy(intersect).^2;
dt(intersect) = dt(intersect) - sqrt(dist2./v2);
part.x(intersect) = ix;
part.y(intersect) = iy;
end

function [part,dt] = checkBoundaries(part,sim,const,dt,forbid,currentBoundary,x,y)
[intersect,ix,iy] = lineIntersection(currentBoundary,[part.x;part.y;x;y],forbid);
if sum(intersect) == 0
    return;
end
[part,dt] = movePointToCrossing(part,dt,intersect,ix(intersect),iy(intersect));
if currentBoundary(1) == currentBoundary(3) %vertical
    if sim.question == 3 && sim.scatter
        part = assignScatterBoundaryVel(part,sim,const,intersect,1);
    else
        part.vx(intersect) = -part.vx(intersect);
    end
elseif currentBoundary(2) == currentBoundary(4) %horizontal
    if sim.question == 3 && sim.scatter
        part = assignScatterBoundaryVel(part,sim,const,intersect,0);
    else
        part.vy(intersect) = -part.vy(intersect);
    end
end
end

function [intersect,ix,iy] = lineIntersection(mainLine,lines,include)
if nargin == 2
    include = true(1,size(lines,2));
end

ix = zeros(1,length(include));
iy = zeros(1,length(include));
intersect = false(1,length(include));

s1x = mainLine(3) - mainLine(1);
s1y = mainLine(4) - mainLine(2);
s2x = lines(3,include) - lines(1,include);
s2y = lines(4,include) - lines(2,include);

s = (-s1y .* (mainLine(1) - lines(1,include)) + s1x .* (mainLine(2) - lines(2,include))) ./ (-s2x .* s1y + s1x .* s2y);
t = ( s2x .* (mainLine(2) - lines(2,include)) - s2y .* (mainLine(1) - lines(1,include))) ./ (-s2x .* s1y + s1x .* s2y);

intersect(include) = s >= 0 & s <= 1 & t >= 0 & t <= 1;
ix(intersect) = mainLine(1) + (t(intersect(include)) * s1x);
iy(intersect) = mainLine(2) + (t(intersect(include)) * s1y);
end

function [part] = assignScatterBoundaryVel(part,sim,const,intersect,orientation)
expectedSigns = zeros(1,sim.particles);
if orientation == 1
    expectedSigns(intersect) = -sign(part.vx(intersect));
elseif orientation == 0
    expectedSigns(intersect) = -sign(part.vy(intersect));
end
IAssign = intersect;
while sum(IAssign) ~= 0
    part = assignVelocities(part,sim,const,IAssign);
    if orientation == 1
        IAssign(intersect) = sign(part.vx(intersect)) ~= expectedSigns(intersect);
    elseif orientation == 0
        IAssign(intersect) = sign(part.vy(intersect)) ~= expectedSigns(intersect);
    end
end
end

function [val] = wrap(val,Min,Max)
val = mod(val,Max - Min) + Min;
end

function [part] = calcValues(part,const,dt)
part.colTime = part.colTime + dt;
part.v2 = part.vx.^2 + part.vy.^2;
part.v = sqrt(part.v2);
part.colDist = part.colDist + dt * part.v;
part.temp = part.v2.*const.meff/(2*const.kb);
end

function [part] = updateAdd(part,lastX,lastY,lengthX)
indexToPlot = part.traceIndex & (abs(lastX - part.x) < (lengthX/2) );
part.addX = [lastX(indexToPlot);part.x(indexToPlot)];
part.addY = [lastY(indexToPlot);part.y(indexToPlot)];
part.addC = part.colour(indexToPlot);
end

function [part] = scatter(part,sim,const)
if sim.question > 1
    I = (rand(1,sim.particles) < sim.Pscat);
    part.colDist(I) = 0;
    part.colTime(I) = 0;
    part = assignVelocities(part,sim,const,I);
end
end

function [stat] = updateStats(stat,part,dt,currentStep)
stat.time(currentStep) = dt*(currentStep-1);
stat.temp(currentStep) = mean(part.temp);
stat.avtemp(currentStep) = mean(stat.temp(1:currentStep));
stat.MFP(currentStep) = mean(part.colDist);
stat.tmn(currentStep) = mean(part.colTime);
stat.avV(currentStep) = mean(part.v);
stat.avVrms(currentStep) = sqrt(mean(part.v2));
[count,temp] = generateTemperatureAverage(stat,part);
stat.avHeatMap = stat.avHeatMap + count;
stat.avTempMap = stat.avTempMap + temp;
end

function [count,temp] = generateTemperatureAverage(stat,part)
count = zeros(length(stat.heatMap.x),length(stat.heatMap.y));
temp = zeros(length(stat.heatMap.x),length(stat.heatMap.y));
for n = 1:(length(stat.heatMap.x)-1)
    for m = 1:(length(stat.heatMap.y)-1)
        I = part.x > stat.heatMap.x(n) & part.x < stat.heatMap.x(n+1) & part.y > stat.heatMap.y(m) & part.y < stat.heatMap.y(m+1);
        count(n,m) = sum(I);
        temp(n,m) = mean(part.temp(I));
    end
end
end

function [plt] = updatePlot(plt,part,stat,sim,currentStep,createPlot)
if nargin == 5
    createPlot = sim.movie;
end

for n = 1:sum(size(part.addX,2))
    line(plt.fig1.ax1,part.addX(:,n)*sim.dist.mult,part.addY(:,n)*sim.dist.mult,'Color',part.addC(n));
end

if ~createPlot
    return;
end

if sim.movie
    set(plt.fig1.plotPos, 'XData', part.x.*sim.dist.mult, 'YData', part.y.*sim.dist.mult);
    set(plt.fig1.plotVel,'XData',part.x.*sim.dist.mult,'YData',part.y.*sim.dist.mult,'UData',part.vx.*sim.dist.mult,'VData',part.vy.*sim.dist.mult);
end

if sim.question ~= 3
    title(plt.fig2.ax1,sprintf('System Temperature (Average = %.0f K)',stat.avtemp(currentStep)));
    set(plt.fig2.plotTemp, 'XData', stat.time(1:currentStep)*sim.time.mult, 'YData', stat.temp(1:currentStep));
    if sim.question ~= 1
        set(plt.fig2.plotAvTemp, 'XData', stat.time(1:currentStep)*sim.time.mult, 'YData', stat.avtemp(1:currentStep)); 
    end
end

if sim.question == 2
    delete(plt.fig3.vHist)
    plt.fig3.vHist = histogram(plt.fig3.ax1,part.v*sim.speed.mult,'FaceColor','g');
    plt.fig3.axis = updateAxis(plt.fig3.axis,axis(plt.fig3.ax1));
    axis(plt.fig3.ax1,plt.fig3.axis);
    set(plt.fig3.Vav,'Value',stat.avV(currentStep)*sim.speed.mult);
    set(plt.fig3.Vrms,'Value',stat.avVrms(currentStep)*sim.speed.mult);
    title(plt.fig3.ax1,sprintf('Particle Velocity Distribution at t = %f %s (Mean = %.0f %s)',stat.time(currentStep)*sim.time.mult,sim.time.unit,stat.avVrms(currentStep)*sim.speed.mult,sim.speed.unit));
end
end

function [out] = updateAxis(oldAxis,newAxis)
out = zeros(1,4);
out(1) = min([oldAxis(1),newAxis(1)]);
out(2) = max([oldAxis(2),newAxis(2)]);
out(3) = min([oldAxis(3),newAxis(3)]);
out(4) = max([oldAxis(4),newAxis(4)]);
end

function [] = outputCalcResults(const,sim,stat)
    fprintf('vth = %.0f km/s\n',const.vth/1000);
    fprintf('vav = %.0f km/s\n',const.vav/1000);
    fprintf('Mean Free Path = %.3f nm\n',const.MFP*1e9);
    fprintf('Total Sim Time = %.0f ps\n',sim.dt*sim.steps*1e12);
    fprintf('Probability of Scatter = %e \n',sim.Pscat);    
    fprintf('Probability of Scatter = %e \n',sim.Pscat);
    fprintf('Mean Free Path (Experimental) = %f nm\n',getPlateauAverage(stat.MFP*1e9));
    fprintf('Mean Time Between Collisions (Experimental) = %f ps\n',getPlateauAverage(stat.tmn*1e12));
    fprintf('Mean Free Path = %.3f nm\n',mean(stat.avV)*const.tmn*1e9);
end

function [out] = getPlateauAverage(y)
y = flip(y);
s = 0;
for n = 1:length(y)
    s = s + y(n);
    if abs(s/n - y(1)) > y(1)/100
        out = s/n;
        fprintf('Average Over %d Points\n',n);
        return;
    end
end
end

function [] = plotHeatMaps(board,sim,stat)
if sim.question == 3
plotParticleHeatMap(board,sim,stat);
plotTempHeatMap(board,sim,stat);
end
end

function [fig] = plotParticleHeatMap(board,sim,stat)
fig.h = figure;
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
axis(fig.ax1,[board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
pcolor(fig.ax1,stat.heatMap.X,stat.heatMap.Y,(stat.avHeatMap')/sim.steps);
colorbar(fig.ax1);
xlabel(fig.ax1,sprintf('x (%s)',sim.dist.unit));
ylabel(fig.ax1,sprintf('y (%s)',sim.dist.unit));
title(fig.ax1,'Particle Distribution HeatMap');
grid(fig.ax1,'on');
hold(fig.ax1, 'on')
end

function [fig] = plotTempHeatMap(board,sim,stat)
fig.h = figure;
fig.ax1 = subplot(1, 1, 1, 'Parent', fig.h);
axis(fig.ax1,[board.xMin,board.xMax,board.yMin,board.yMax]*sim.dist.mult);
pcolor(fig.ax1,stat.heatMap.X,stat.heatMap.Y,(stat.avTempMap')/sim.steps);
colorbar(fig.ax1);
xlabel(fig.ax1,sprintf('x (%s)',sim.dist.unit));
ylabel(fig.ax1,sprintf('y (%s)',sim.dist.unit));
title(fig.ax1,'Temperature Distribution HeatMap');
grid(fig.ax1,'on');
hold(fig.ax1, 'on')
end
