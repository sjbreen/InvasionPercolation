%This is the control script for the 2D invasion percolation routine

clear all; clc

vidnovid = 0; %0 for no video, 1 for video
barriernobarrier = 0; %0 for no barrier, 1 for mid chamber barrier

Rmax = 0.0001; %maximum pore radius [m]
Rmin = 0.00001; %minimum pore radius [m]
FineRmax = 0.00001;
FineRmin = 0.000001;
IntTen = 0.05; %interfacial tension -> 0.05 [N/m] for air/water
ConAng = 110; %contact angle [deg] measured through the invading phase. B/t 90 and 180 if non-wetting
DensDef = 1000; %density of defending fluid [kg/m^3]
DensInv = 1.2; %density of invading fluid [kg/m^3]
Grav = 9.8; %acceleration of gravity [m/s^2]
Height = 0.57; %height of medium [m]
Width = 0.28; %width of medium [m]
fractionclosed = 0.5; %fraction of radii equal to zero

%Define a capillary pressure field. The buoancy forces need to be
%calculated dynamically, during the local invasion "search".

RField = abs(normrnd(((Rmax+Rmin)/2),((Rmax-Rmin)/4),Height*1e3,Width*1e3)); %build a normal random field of radii
maxind = RField>=Rmax;
RField(maxind) = Rmax;

if barriernobarrier == 1;
    RField(1:size(RField,1)/2,:) = abs(normrnd(((FineRmax+FineRmin)/2),((FineRmax-FineRmin)/4),Height*1e3/2,Width*1e3));
end

%Make a random collection of radii equal to zero (note that there may be
%some repeating i,j pairs, making actual fractionclosed slightly smaller
numclosed = fractionclosed*numel(RField);
randi = round(1+(size(RField,1)-1)*rand(numclosed,1));
randj = round(1+(size(RField,2)-1)*rand(numclosed,1));
for i = 1:size(randi,1);
    RField(randi(i),randj(i)) = 0;
end

CapPres = -2*IntTen*cosd(ConAng)*(1./RField); %capillary pressure field
A = CapPres;

%Run Invasion Percolation routine and imshow results.

StartPoint = [floor(0.95*size(RField,1)),floor(0.5*size(RField,2))]; %Gas
% injection point in A, 90% down, 50% across
BW = 5; %Border Width

if vidnovid == 1
    [A,InvNum,InvList,aviobj] = InvPFuncListed(A,StartPoint,BW,DensDef,DensInv,Grav,CapPres,vidnovid);
elseif vidnovid == 0
    [A,InvNum,InvList] = InvPFuncListed(A,StartPoint,BW,DensDef,DensInv,Grav,CapPres,vidnovid);
    PathSub = A-CapPres; %This will make blocked paths = NaN and fluid path = InvNum (=Inf)
    Path = PathSub>=InvNum; %Pathway is = Inf
    imshow(Path);
end

%Correct #s for experiments: Rmax: 0.0001, Rmin: 0.00001, IntTen = 0.05,
%ConAng = 110, DensDef = 1000, DensInv = 1.2, Grav = 9.8, Height = 0.57,
%Width = 0.28, fractionclosed = 0.5