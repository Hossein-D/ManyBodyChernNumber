# A matlab code to evalulate the Chern number using the swap formula introduced in PhysRevB.103.075102 for Chern insulators. 
function [numeratorSwap, denominatorSwap] = MCChernInsCylinder(state, Nx, Ny, NMCS, plotOn)

% NMCS, theta, swapLen

% Cylinder
% Y-Direction Periodic 

NTheta = 100;
theta = 0:1/NTheta:2*pi;
Ns = Nx * Ny;
Np = 1/2 * Ns;
latticeConst = 1;

Lx = Nx * latticeConst;
Ly = Ny * latticeConst;

% Monte Carlo Hopping Distance
MCDelta = 1;

grid = zeros(2, Ns);
grid(1, 1:Ns) = floor((0:Ns-1)/Ny);
grid(2, 1:Ns) = (0:Ns-1) - Ny*floor((0:Ns-1)/Ny);

vectorGrid = 0:Ns-1;

if nargin < 5
    theta  = pi;
    swapLen = 1;
end

% Initialization
OccPosit = zeros(2, Np);
initOccIndex = randsample(Ns, Np)-1;
OccIndex = zeros(Ns, 1);

for n = 1:Ns
    if(find(initOccIndex==n-1))
        OccIndex(n) = 1;
    end
end
OccIndex';

initOccPosit(1, 1:Np) = floor(initOccIndex(1:Np)/Ny);
initOccPosit(2, 1:Np) = initOccIndex(1:Np) - Ny * floor(initOccIndex(1:Np)/Ny);
OccPosit = initOccPosit;
    
initOccPosit;

if plotOn ==1 
    figure()
    plot(initOccPosit(1, :), initOccPosit(2, :), '.k', 'MarkerSize',30)    
end
removedOccPosit = zeros(2, Np-1);

randomProb = rand(NMCS*Np, 1);
    
randomGenAngleVec = randi(4, [NMCS*Np, 1]) - 1;
oldSlaterMatrix = zeros(Np, Np);
newSlaterMatrix = zeros(Np, Np);

oldSiteIndexVec = zeros(1, Np);
newSiteIndexVec = zeros(1, Np);
%
for MCS = 1:NMCS    
    acceptance = 0.0;
    for n = 1 : Np
        
        oldOccX = OccPosit(1, n);
        oldOccY = OccPosit(2, n);
        
        
        randAngle = pi/2*randomGenAngleVec(Np*(MCS-1) + n);
        newOccX = round(oldOccX + MCDelta * cos(randAngle), 0);
        newOccY = round(mod(oldOccY + MCDelta * sin(randAngle), Ly), 0);
        %[newOccX, newOccY]
        
        oldSiteIndex = oldOccY + Ny * oldOccX + 1;
        oldSiteIndex = round(oldSiteIndex, 0);
        
        newSiteIndex = round(newOccY + Ny * newOccX + 1, 0);
        
        if (newOccX < 0 || newOccX >= Lx)
            %display('continue')
            continue
        end
        
        if OccIndex(newSiteIndex)
            continue
        end
        oldSiteIndexVec(1, 1:Np) = round(OccPosit(2, 1:Np) + Ny * OccPosit(1, 1:Np) + 1, 0);
        newSiteIndexVec(1, 1:n-1) = OccPosit(2, 1:n-1) + Ny * OccPosit(1, 1:n-1) + 1;
        newSiteIndexVec(1, n+1:Np) = OccPosit(2, n+1:Np) + Ny * OccPosit(1, n+1:Np) + 1;
        
        newSiteIndexVec(1, n) = newSiteIndex;
        
        newSiteIndexVec = round(newSiteIndexVec, 0);
        %newSiteIndex = newOccY + Ny * newOccX + 1;
        %newSiteIndex = round(newSiteIndex, 0);
        
        %removedOccPosit(1, :) = [OccPosit(1, 1:n-1) OccPosit(1, n+1:Np)];
        %removedOccPosit(2, :) = [OccPosit(2, 1:n-1) OccPosit(2, n+1:Np)];        
        %oldSiteIndexVec;
        oldSlaterMatrix(:, 1:Np) = state(oldSiteIndexVec, 1:Np);
        %newSiteIndexVec;
        newSlaterMatrix(:, 1:Np) = state(newSiteIndexVec, 1:Np);
        
        oldProb = abs(det(oldSlaterMatrix)/sqrt(factorial(Np)));
        newProb = abs(det(newSlaterMatrix)/sqrt(factorial(Np)));
        
        %logDistanceNew = sum(log(abs(exp(1i*gamma*newCompOccPosit) - exp(1i*gamma*removedOccPosit(:))).^2));        

%        logDistanceNew = sum(log((newWX - Position(1,0:n)).^2 + ...
%                       (newWY - removedPosition(2,0:n)).^2) + ...
%                   (newWX - Position(1,n+1:N)).^2 + ...
%                       (newWY - removedPosition(2,n+1:N)).^2) ...
        %print(logDistanceNew)
        probRatio = newProb/oldProb;
        
        %newSiteIndex;
        if probRatio > randomProb(Np*(MCS-1) + n)
            %display('move');
            OccIndex(oldSiteIndex) = 0;
            OccIndex(newSiteIndex) = 1;
            OccPosit(1, n) = newOccX;
            OccPosit(2, n) = newOccY;
        end
    end
end


if plotOn == 1
    hold on
    plot(OccPosit(1, :), OccPosit(2, :), 'or', 'MarkerSize',15) 
end

ell0 = 1;
swapLen = 3*ell0;

if Lx >= 10*ell0
    L0 = 3*ell0;
    L1 = swapLen + L0;
    L2 = swapLen + L1;
    L3 = swapLen + L2;
else
    return;
end


swapOccPosit = zeros(2, Np);
swapOccPosit(2, :) = OccPosit(2, :);


Wtau = 1;
Wchi = ones(1, NTheta);

for n = 1 : Np    
    if OccPosit(1, n) > L0 && OccPosit(1, n) <= L2
        Wtau = Wtau * exp(2*pi*1i*OccPosit(2, n)/Ly);
    elseif OccPosit(1, n) > L1 && OccPosit(1, n) <= L3
        Wchi(1:NTheta) = Wchi(1:NTheta) .* exp(1i * theta(1:NTheta));
    end
    
    if OccPosit(1, n) > L0 && OccPosit(1, n) <= L1
        swapOccPosit(1, n) = OccPosit(1, n) + 2*swapLen;
    elseif OccPosit(1, n) > L2 && OccPosit(1, n) <= L3
        swapOccPosit(1, n) = OccPosit(1, n) - 2*swapLen;
    else 
        swapOccPosit(1, n) = OccPosit(1, n);
    end
end


if plotOn==1
    plot(swapOccPosit(1, :), swapOccPosit(2, :), '.g', 'MarkerSize',30) 
    th = 0:pi/50:2*pi;
    xunit0 = th - th + L0;
    yunit0 = th;
    plot(xunit0, yunit0, 'Markersize', 4);

    xunit1 = th - th + L1;
    yunit1 = th;
    plot(xunit1, yunit1, 'Markersize', 4);

    xunit2 = th - th + L2;
    yunit2 = th;
    plot(xunit2, yunit2, 'Markersize', 4);

    xunit3 = th - th + L3;
    yunit3 = th;
    plot(xunit3, yunit3, 'Markersize', 4);
end        

oldSiteIndexVec(1, 1:Np) = round(OccPosit(2, 1:Np) + Ny * OccPosit(1, 1:Np) + 1, 0);
swapSiteIndexVec(1, 1:Np) = swapOccPosit(2, 1:Np) + Ny * swapOccPosit(1, 1:Np) + 1;

SlaterMatrix(:, 1:Np) = state(oldSiteIndexVec, 1:Np);
swapSlaterMatrix(:, 1:Np) = state(swapSiteIndexVec, 1:Np);

PsiStar = conj(det(SlaterMatrix)/sqrt(factorial(Np)));
swapPsiStar = conj(det(swapSlaterMatrix)/sqrt(factorial(Np)));

swapPsiOverPsi_Star = swapPsiStar/PsiStar;

numeratorSwap = zeros(1, NTheta);
denominatorSwap = zeros(1, NTheta);

numeratorSwap = Wtau*Wchi*swapPsiOverPsi_Star;
denominatorSwap = Wchi*swapPsiOverPsi_Star;

end
