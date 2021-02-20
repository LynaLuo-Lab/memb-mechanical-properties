%-------------------------------------------------------------------------
% LTF Method for leaflet Ka moduli
%
% Method based on Doktorova et al. 2019, Biophysical Journal 116, 487-502
%
% Milka Doktorova, March 2019
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% calculate_Ka.m
%
% Parses the heightmap.dat and atomGroups.txt files to create necessary
% input parameters for leaflet Ka calculation.
%
% Calculates error on the leaflet Ka with 2D-bootstrapping.
% NOTE: uses parallel computing toolbox in MATLAB. If not available,
% replace 'parfor' with 'for'.
%-------------------------------------------------------------------------

%-------------------------------------------------------------------------
% Input parameters
%
% leaf      1 for top, 0 for bottom
% frame     frame # from which the analysis should start
%           (e.g. use when you want to discard the first N frames)
% boot_flag calculate errors: 1 for yes, 0 for no
% db1       carbon position/s of double bond/s in first chain (i.e. chain
%           that appears right after headgroup selection in atomGroups.txt)
%           e.g. 0 for DMPC, 9 for POPC, [5 8 11 14] for SAPE
% db2       same as in db1 but for second chain (9 for DOPC, 0 for DMPC)
% db3       same as in db1 but for third chain (9 for TOCL, 0 for DMPC)
% db4       same as in db1 but for fourth chain (9 for TOCL, 0 for DMPC)
%
% Example calculation for top leaflet of POPC Ka (in example/ directory):
% ka = calculate_Ka(1,0,1,9,0,0,0)
% -------------------------------------------------------------------------
% Input files
%
% atomGroups.txt
%           same as atomNames.txt file but with empty lines separating
%           the different sets of atoms (headgroup, chain 1, chain 2, etc.)
% apl.txt
%           average area per lipid in the leaflet
% temperature.txt
%           temperature of simulation in Kelvin
% -------------------------------------------------------------------------
% Output files
% 
% avgsn_analysis_xxx
%           intermediate analysis plots including distribution of relevant
%           thicknesses, PMF and fit, and slopes analysis for identifying
%           the relevant surface
%
% Ka_avgsn_xxx
%           correlation vs distance plot, relevant surface and resulting Ka
%
% Ka_avgsn_xxx.txt
%           Ka, carbon at relevant surface, average thickness at relevant
%           surface, range for PMF fit
%
% bootstrap_100_xxx
%           bootstrap information:
%           [{number of bootstraps} {number of frames in a block}
%            {number of grid points in a block} {number of blocks} 
%            {length of re-generated dataset as a fraction of length of
%            original dataset)]
%
% Ka_bootstrap_100_xxx
%           all bootN Ka values calculated from bootstrapping
%
% Ka_stat_100_xxx
%           mean Ka, std and 95% confidence interval
% 
%-------------------------------------------------------------------------
function [mKa] = calculate_Ka(leaf, frame, boot_flag, db1, db2, db3, db4)

%-------------------------------------------------------------------------
% USER SPECIFIED PARAMETERS
%-------------------------------------------------------------------------
% specify the number of the atom in the average chain that forms the 
% reference surface (usually, the first atom not connected to oxygen)
% default: 2 (i.e. the average of the z-pos of the second carbon in both 
% chains, assuming that in atomGroups.txt the chains start with the first 
% carbon on the chain)
headCarbon = 2;

bootN = 100; % number of bootstraps

%-------------------------------------------------------------------------
% read files
a0=load('apl.txt');
temper = load('temperature.txt');
h = load('heightmap.dat');

% read file with atom names where common (i.e. headgroup/backbone) atoms 
% and chains are separated by empty line
at = textread('atomGroups.txt','%s','delimiter','\r\n');
% find empty lines in file that separate common atoms and chains
ind = strmatch(' ',at);

elements = zeros(length(ind)+1,1);

elements(1) = ind(1)-1;

for i=1:length(ind)-1
    elements(i+1) = length(ind(i)+1:ind(i+1)-1);
end
elements(end) = length(ind(end)+1:length(at));

h = h(h(:,1)>=frame,:);

% make z-positions positive if bottom leaflet
if leaf==0
    h = -h;
    leaflb = 'bot';
else
    leaflb = 'top';
end

% identify the size of the grid map and number of frames
i = 1;
while h(i,2) == h(1,2)
    i = i+1;
end
ngridX = i-1;
ngridY = length(h(h(:,1)==h(1,1)))/ngridX;
ngrids = ngridX*ngridY;
nfr = length(h(:,1))/ngrids;

% trim the frame and grid info, and the common atoms
h = h(:,4+elements(1):end);
nchains = length(elements)-1;

% correct for presence of double bonds
% see paragraph after Eq. 15 in paper
dbflag = zeros(1,length(h(1,:)));
if db1>0
    dbflag(db1) = 1;
end
if db2>0
    dbflag(elements(2)+db2) = 1;
end
if nchains>2
    dbflag(elements(2)+elements(3)+db3) = 1;
    dbflag(elements(2)+elements(3)+elements(4)+db4) = 1;
end

hh = [];
i = 1;
while i<=length(h(1,:))
    if i<length(h(1,:)) && dbflag(i+1)==1
        hh = [hh (1/3)*sum(h(:,i:i+2),2)];
        i = i+3;
    else
        hh = [hh h(:,i)];
        i = i+1;
    end
end
elements(2) = elements(2)-length(db1(db1>0))*2;
elements(3) = elements(3)-length(db2(db2>0))*2;
if nchains>2
    elements(4) = elements(4)-length(db3(db3>0))*2;
    elements(5) = elements(5)-length(db4(db4>0))*2;
end
h = hh;

% calculate the z-pos of an average chain
nDataPoints = length(h(:,1));
maxChainLength = max(elements(2:end));
sn = zeros(nDataPoints,maxChainLength);
startInd = 1;

for i=1:nchains
    curr_chain = h(:,startInd:startInd+elements(i+1)-1);
    curr_chainLength = length(curr_chain(1,:));
        
    if maxChainLength > curr_chainLength
        sn = sn + (1/nchains)*[curr_chain(:,1:curr_chainLength) curr_chain(:,end)*ones(1,maxChainLength-curr_chainLength)];
    else
        sn = sn + (1/nchains)*curr_chain;
    end
    
    startInd = startInd + elements(i+1);
end

mKa = Ka_(sn(:,headCarbon:end),a0,temper,1,leaflb,headCarbon);


%-------------------------------------------------------------------------
% ERROR CALCULATION WITH 2D BOOTSTRAP
%-------------------------------------------------------------------------
if ~boot_flag
    return
end

front_thick = sn(:,headCarbon:headCarbon+2)-sn(:,end);
frontT = length(front_thick(1,:))


% Find optimal time segment
% Consider only first 3 atoms since they'll likely have the slowest
% autocorrelation times
parfor iddif=1:frontT
    th = front_thick(:,iddif);
    disp(iddif);
    % find stride s.t. abs(acf) < 0.2
    stride = 1;
    ss = zeros(ngrids,1);
    for pts=1:ngrids
        [y rv] = acf(th(pts:ngrids:end),floor(nfr/2));
        ff = find(y<=rv);
        if ff(1)>stride
            stride = ff(1);
        end       
        ss(pts) = ff(1);
    end
    dlmwrite(strcat('stride_',num2str(iddif),'.dat'),ceil(median(ss)),'delimiter',' ');
    disp(ceil(median(ss)));
end

timeN = 0;
for iddif = 1:frontT
    tN = load(strcat('stride_',num2str(iddif),'.dat'));
    tN = tN(1);
    if tN > timeN
        timeN = tN;
    end
end

% Calculate error

spaceN = floor(min(ngridX,ngridY)/2)^2;
ntimes = floor(nDataPoints/(spaceN*timeN));
chunk = spaceN*timeN;

Ka_data = [];
parfor boot = 1:bootN
    newsn = zeros(timeN*spaceN*ntimes,1);
    for i=1:ntimes
        stFr = randi(nfr-timeN+1);
        frames = stFr:stFr+timeN-1;
        cornerPtx = randi(ngridX) + ngridX;
        cornerPty = randi(ngridY) + ngridY;
        dummypoints = reshape(1:ngrids,[ngridY,ngridX]);
        dummyrep = [dummypoints dummypoints dummypoints;dummypoints dummypoints dummypoints;dummypoints dummypoints dummypoints];
        nghpts = reshape(dummyrep(cornerPty:cornerPty+sqrt(spaceN)-1,cornerPtx:cornerPtx+sqrt(spaceN)-1),[spaceN 1]);
        nghpts = sort(nghpts);
        sntemp = [];
        for j=1:spaceN
            sntemp = [sntemp; [(stFr-1)*ngrids+nghpts(j):ngrids:(stFr+timeN-2)*ngrids+nghpts(j)]'];
        end
        newsn((i-1)*chunk+1:i*chunk) = sntemp;
    end
    
    newmKa = Ka_(sn(newsn,headCarbon:end),a0,temper,0,leaflb,headCarbon);

    disp([boot newmKa])
    Ka_data = [Ka_data; newmKa];
end

dlmwrite(strcat('Ka_bootstrap_',num2str(bootN),'_',leaflb,'.dat'),Ka_data,' ');
dlmwrite(strcat('bootstrap_',num2str(bootN),'_',leaflb,'.dat'),[bootN timeN spaceN ntimes 100*(chunk*ntimes)/(nfr*ngrids)],' ');
Kasort = sort(Ka_data);
data = [mKa std(Ka_data) Kasort(0.05*bootN) Kasort(0.95*bootN)];
dlmwrite(strcat('Ka_stat_',num2str(bootN),'_',leaflb,'.txt'),data,' ');
disp([leaflb ' Ka: ' num2str(mKa) ' +/- ' num2str(std(Ka_data))])

return

