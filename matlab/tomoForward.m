function Cflat = tomoForward(mode, step, ph, numDist)
%TOMOFORWARD Simulate concentration and column density profiles of plumes
%   Jan 2024
close all;

%wasm.console.log('Starting tomoForward');
%wasm.console.log(sprintf('Mode is %f.', mode));
%wasm.console.log(sprintf('Step is %f.', step));
%wasm.console.log(sprintf('Plume height is %f.', ph));
%wasm.console.log(sprintf('numDist is %f.', numDist));

assert(isa(mode, 'double'));
assert(all( size(mode) == [ 1 ]));

assert(isa(step, 'double'));
assert(all( size(step) == [ 1 ]));

assert(isa(ph, 'double'));
assert(all( size(ph) == [ 1 ]));

assert(isa(numDist, 'double'));
assert(all( size(numDist) == [ 1 ]));

% Determine the scanning mode
%mode = 1; % scan: 1, traverse: 2
%step = 3.6; % scanning step for scan mode
%ph = 1; % plume height scale
%numDist = 3; % 0 for one distribution, >0 for mixed distribution

% Compute the pdf using a combination of bivariate normal distributions
% Create a grid of evenly spaced points in two-dimensional space
x = -5*ph:0.1*ph:5*ph;
y = 0*ph:0.1*ph:3*ph;
[X,Y] = meshgrid(x,y);
P = [X(:) Y(:)];

% Define the parameters of the distributions
mu = cell(1,numDist);
Sigma = cell(1,numDist);
corr = zeros(1,numDist);
C0 = mvnpdf(P,[randi([-1 1],1) ph],[0.1 0;0 0.1]);
C0 = reshape(C0,length(y),length(x));
for i = 1:numDist
    mu{i} = [2*ph*rand(1,1) ph*rand(1,1)+1];
    corr(i) = randi([-1,1])*rand(1,1)/10;
    Sigma{i} = [rand(1,1)+0.1 corr(i); corr(i) rand(1,1)+0.01];
end
for i = 1:numDist
    C1 = mvnpdf(P,[mu{i}],[Sigma{i}]);
    C1 = reshape(C1,length(y),length(x));
    C0 = C0+C1;
end
for i = 1:length(y)
    for j = 1:length(x)
        if C0(i,j) < 1e-1*max(max(C0))
            C0(i,j) = 0;
        end
    end
end
C = C0/sum(sum(C0));

% Compute the column density along the viewing path and determine plume
% parameters
coder.varsize('alpha');
alpha = -90:step:90;
if mode>1
    alpha = zeros(1,length(alpha));
end
Cs = 0*C;
coder.varsize('C')
S = zeros(length(alpha),1);
V = S;
if mode<2
    for m = 1:length(alpha)
        for i = 1:length(x)
            for j = 1:length(y)
                if abs(atand(x(i)/y(j))-alpha(m))<=mean(diff(alpha))
                    Cs(j,i) = C(j,i);
                    S(m) = S(m)+Cs(j,i);
                    V(m) = S(m)*cosd(alpha(m));
                    Cs(j,i) = 0;
                end
            end
        end
    end
    deltaX = abs(diff(tand(alpha)));
    Vm = V(1:end-1) + diff(V)/2;
    deltaS = Vm(2:end-1).*deltaX(2:end-1)';
    plumeCentre = dot(V,alpha)/sum(V);
    plumeCompleteness = 1 - 0.5*(max(0.1*sum(V(1:5)),0.1*sum(V(end-4:end)))/max(S));
    [valueMax,posMax] = max(Vm);
    posWidth = find((Vm/valueMax)<(0.1));
    wLeft = max(find((posWidth-posMax)<0));
    valueLeft = posWidth(wLeft);
    wRight = min(find((posWidth-posMax)>0));
    valueRight = posWidth(wRight);
    plumeWidth = abs(tand(alpha(valueRight))-tand(alpha(valueLeft)));
else
    for i = 1:length(x)
        S(i) = sum(sum(C(:,i)));
    end
    V = S;
    deltaX = abs(diff(x));
    Vm = V(1:end-1) + diff(V)/2;
    deltaS = Vm(2:end-1).*deltaX(2:end-1)';
    plumeCentre = dot(V,x)/sum(V);
    plumeCompleteness = 1 - 0.5*(max(0.1*sum(V(1:5)),0.1*sum(V(end-4:end)))/max(S));
    [valueMax,posMax] = max(Vm);
    posWidth = find((Vm/valueMax)<(0.1));
    wLeft = max(find((posWidth-posMax)<0));
    valueLeft = posWidth(wLeft);
    wRight = min(find((posWidth-posMax)>0));
    valueRight = posWidth(wRight);
    plumeWidth = abs(x(valueRight)-x(valueLeft));
end

% Plot the scanning pattern
subplot(211);
beta = 90-alpha;
for m = 1:length(alpha)
    line([0 max(max(x),max(y))*cosd(beta(m))],[0 max(max(x),max(y))*sind(beta(m))],'linestyle','-','color',[0.5 0.5 0.5]);
end
hold on;

% Plot the concentration profile
contour(x,y,C);
xlabel('x');
ylabel('y');
colorbar('southoutside');
grid on;
title('Concentration');

% Plot the column densities
subplot(212)
if mode<2
    plot(alpha,S,'-ok',alpha,V,'r');
    hold on;
    line([plumeCentre plumeCentre],[0 1.5*valueMax],'linestyle','--','color','b');
    line([alpha(valueLeft) alpha(valueLeft)],[0 1.5*valueMax],'linestyle','--','color',[0.5 0.5 0.5]);
    line([alpha(valueRight) alpha(valueRight)],[0 1.5*valueMax],'linestyle','--','color',[0.5 0.5 0.5]);
    xlabel('view angle / deg');
    xlim([-90 90]);
else
    plot(x,S,'-ok',x,V,'r');
    hold on;
    line([plumeCentre plumeCentre],[0 1.5*valueMax],'linestyle','--','color','b');
    line([x(valueLeft) x(valueLeft)],[0 1.5*valueMax],'linestyle','--','color',[0.5 0.5 0.5]);
    line([x(valueRight) x(valueRight)],[0 1.5*valueMax],'linestyle','--','color',[0.5 0.5 0.5]);
    xlabel('x');
end
legend('slant','vertical','centre of mass','limits');
ylabel('column density');
%title(['Max pos: ',num2str(plumeCentre,'%.1f'),' [deg or heights]; Completeness: ',num2str(plumeCompleteness,'%.1f'),'; Width: ',num2str(plumeWidth,'%.1f'),' [heights]']);
grid on;

Cflat = reshape(C,1,3131);

clc;
end

