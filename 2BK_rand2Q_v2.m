options = optimset('Display','off');
%% Load Precipitation Time-series %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for RainTSinputSelection = 1:3
if RainTSinputSelection == 1
%% Time-Series #1
    load('SmithRivP.mat'); %88-02 (15 years - Leap yrs 1,5,9,13)
    P_t = SmithRivP(14611:20089); %Precip TS + experiment
    
    leapYearGuide=[365;365;365;366;365;365;365;366;365;365];
    doy = [1:1:leapYearGuide(1),1:1:leapYearGuide(2),1:1:leapYearGuide(3),...
        1:1:leapYearGuide(4),1:1:leapYearGuide(5),1:1:leapYearGuide(6),...
        1:1:leapYearGuide(7),1:1:leapYearGuide(8),1:1:leapYearGuide(9),1:1:leapYearGuide(10)];
    TScutoff = 1828; %amount of time-series to cut off for spin up (5 years)
    figNums = [1 2];
elseif RainTSinputSelection == 2
%% Time-Series #2
    load('BroadRivP.mat'); %88-02 (15 years - Leap yrs 1,5,9,13)
    P_t = BroadRivP(14611:20089); %Precip TS + experiment
    
    leapYearGuide=[365;365;365;366;365;365;365;366;365;365];
	doy = [1:1:leapYearGuide(1),1:1:leapYearGuide(2),1:1:leapYearGuide(3),...
        1:1:leapYearGuide(4),1:1:leapYearGuide(5),1:1:leapYearGuide(6),...
        1:1:leapYearGuide(7),1:1:leapYearGuide(8),1:1:leapYearGuide(9),1:1:leapYearGuide(10)];
    TScutoff = 1828; %amount of time-series to cut off for spin up (5 years)
    figNums = [3 4];
elseif RainTSinputSelection == 3
%% Time-Series #3
    load('BisleyRivP.mat'); %94,94-06,09 (15 years - Leap yrs 4,8,12)
    P_t = BisleyP; %Precip TS + experiment
    
	leapYearGuide=[365;365;366;365;365;365;366;365;365;365];
	doy = [1:1:leapYearGuide(1),1:1:leapYearGuide(2),1:1:leapYearGuide(3),...
        1:1:leapYearGuide(4),1:1:leapYearGuide(5),1:1:leapYearGuide(6),...
        1:1:leapYearGuide(7),1:1:leapYearGuide(8),1:1:leapYearGuide(9),1:1:leapYearGuide(10)];
    TScutoff = 1827; %amount of time-series to cut off for spin up (5 years)
    figNums = [5 6];
  
else
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hydrologic Model Run
Pavg = mean(P_t(1:end)); %average daily precipitation for Precip TS
ntri = 5000; %number of montecarlo trials
hydroSummary = zeros(ntri,7311); %pre-allocation
parfor tri = 1:ntri %% PARALLEL
    %% KBM STORAGE - FLUX parameter choices
    eta = random('uniform',0.01,0.9);
    SuRef = 10.^(random('Uniform',log10(20),log10(500)));
    SlRef = 10.^(random('Uniform',log10(500),log10(10000)));
    bu = random('uniform',1,20);
    ku = Pavg/(SuRef^bu);
    bl = random('uniform',1,50);
    kl = ((1-eta)*Pavg)/(SlRef^bl);
    dt = 1; %timestep
    tlength = length(P_t); %simulation length (days)
    % Conservative Tracer Time Series
    amp = -10; %seasonal tracer amplitude (permil)
    phs = -4.4; %tracer phase lag (days) 
    %set to be in phase with Smith, which means min Fyw
    avg = 0; %mean tracer value (permil)
    rainTracerFunc = @(t) amp*sin((t*2*pi/365)-phs)+avg;
    Cp_t = rainTracerFunc(1:1:tlength)'; %tracer time-series
    %% KBM STORAGE - FLUX model initialization
    % Pre-allocation & initialization of storage & discharge vectors
    Sl_t = ones(tlength,1); %lower box storage
    Su_t = ones(tlength,1); %upper box storage
    Su_t(1) = SuRef; %initial set to reference storage upper box
    Sl_t(1) = SlRef; %initial set to reference storage lower box
    Ql_t = ones(tlength,2); %lower box discharge
    Qu_t = ones(tlength,3); %upper box discharge (total)
    Lr_t = ones(tlength,1); %upper box discharge (to lower box)
    Cu_t = ones(tlength,2); %upper box tracer concentrations (box, outflow)
    Cu_t(1,2) = avg; Cu_t(1,1) = avg; % initialization 
    Cl_t = ones(tlength,2); %lower box tracer concentration
    Cl_t(1) = avg; %(initialization)
    %% KBM STORAGE - FLUX model Solution Scheme
        for m = 2:length(P_t) % for each day of precipitation    
    % Previous time step values
            Sl = Sl_t(m-1); %lower storage
            Su = Su_t(m-1); %upper storage
            P =  P_t(m); % should this be (previous) or current precipitation?
            %P =  P_t(m-1); %should this be previous or (current)?
            Cu = Cu_t(m-1,1); %upper tracer concentration
            Cl = Cl_t(m-1,1); %lower tracer concentration
	% Upper Box Solution
	% stability check
            stcU = 0.5 + 0.5*(((P-ku*(Su^bu))*dt)/( (P/ku)^(1/bu) - Su));
            pu = min([stcU,1]);%stability criteria
            upperBox = @(xu) (dt * (P - pu*ku*(xu.^bu) - ((1-pu)*ku*(Su.^bu) ))) - (xu - Su);
            Su_t(m) = fsolve(upperBox,0,options); %upper box storage
            L = P + ((Su - Su_t(m))/dt); %upper box outflow
            Lr = (1-eta)*L; %upper box outflow to lower box
            Lr_t(m) = Lr; %save value
            Qu_t(m,1) = eta*L; %upper box outflow to discharge
            Qu_t(m,2) = L; %upper box TOTAL outflow
            Cu_t(m,1) = Cp_t(m,1) + (Cu-Cp_t(m,1)) .* (Su./Su_t(m)).^(P./(P-L)); %[tracer] UPPER  box
            Cu_t(m,2) = ( (P.*Cp_t(m,1)) + (Cu.*Su) - (Cu_t(m,1).*Su_t(m)) )./L ; %[tracer] UPPER  outflow
    % Lower Box
	% stability check
            stcL = 0.5 + 0.5*(((Lr-kl*(Sl^bl))*dt)/( (Lr/kl)^(1/bl) - Sl));
            pl = min([stcL,1]);%stability criteria
            lowerBox = @(xl) (dt * ( Lr - pl*kl*(xl.^bl) - (1-pl)*kl*(Sl.^bl) )) - (xl - Sl);
            Sl_t(m) = fsolve(lowerBox,0,options); %lower box storage
            Ql_t(m,1) = Lr + ((Sl - Sl_t(m))/dt); %lower box outflow
            Ql_t(m,2) = Qu_t(m,1) +  Ql_t(m,1); %TOTAL discharge
            Cl_t(m,1) = Cu_t(m,2) + (Cl-Cu_t(m,2)) .* (Sl./Sl_t(m)).^(Lr./(Lr-Ql_t(m,1))); %[tracer] LOWER box
            Cl_t(m,2) = ( (Lr.*Cu_t(m,2)) + (Cl.*Sl) - (Cl_t(m,1).*Sl_t(m)) )./Ql_t(m,1) ; %[tracer] UPPER  outflow
        end
    %% Model Output Summary
    % Calculate Tracer Time-series from mass balance
    isoQ = ( (Qu_t(:,1).*Cu_t(:,2)) + (Ql_t(:,1).*Cl_t(:,2)) )./ Ql_t(:,2);
    % Q Data saving
    inputParameters = [eta,SuRef,SlRef,bu,ku,bl,kl];
    modelOutput = [Ql_t(TScutoff:end,2)',isoQ(TScutoff:end)'];
    hydroSummary(tri,:) = [inputParameters,modelOutput];
end
%% Extract Working Dataset
%remove param values and clip model spin up period
Qdata =  hydroSummary(:,8:3659); %discharge data
isodata = hydroSummary(:,3660:end); %isotracer data
%zscore data
isodataZ = zscoreWholeMatrix(isodata);
QdataZ =  zscoreWholeMatrix(Qdata);

QdayMean= zeros(length(Qdata(:,1)),366); %pre-allocation
isodayMean= zeros(length(isodata(:,1)),366);  %pre-allocation

for m = 1:366 %for each day
    dayInd = find(doy==m); %find where data from that day is
    QdayMean(:,m) = mean(Qdata(:,dayInd),2); %average
    isodayMean(:,m) = mean(isodata(:,dayInd),2); %average
end
isodataZmean = zscoreWholeMatrix(isodayMean);
QdataZmean =  zscoreWholeMatrix(QdayMean);

%% k-means clustering
nClust = 10;
[IDX,C] = kmeans([QdataZmean,isodataZmean],nClust); %Cluster

%% Find Parameters closest to centroids
clusteredParams = zeros(nClust,7);
for m = 1:nClust
    clustInd = find(IDX==m);
    SSEwiClust= sum(([QdataZmean(clustInd,:),isodataZmean(clustInd,:)] - C(m,:)).^2,2);
    minInd = find(SSEwiClust == min(SSEwiClust));
    clusteredParams(m,:) = hydroSummary(clustInd(minInd),1:7);
end
save(strcat('clusteredQParams',num2str(RainTSinputSelection),'.mat'),...
    'clusteredParams')
clearvars -except options 
end