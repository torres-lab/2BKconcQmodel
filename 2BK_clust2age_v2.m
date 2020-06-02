options = optimset('Display','off');
ndec = 3; %numbder of decimal places to round to

%% Load Clustered Paramater Values & Precipitation Time-series
load('clusteredQParams1.mat')
SRp = clusteredParams;
load('clusteredQParams2.mat')
BRp = clusteredParams;
load('clusteredQParams3.mat')
PRp = clusteredParams;

HySumclust = [SRp;BRp;PRp];
clearvars -except HySumclust options ndec

load('SmithRivP.mat'); %88-02 (15 years - Leap yrs 1,5,9,13)
load('BroadRivP.mat'); %88-02 (15 years - Leap yrs 1,5,9,13)
load('BisleyRivP.mat'); %94,94-06,09 (15 years - Leap yrs 4,8,12)

rainTS(:,1) = SmithRivP(14611:20089); %Precip TS
rainTS(:,2) = BroadRivP(14611:20089); 
rainTS(:,3) = [0;BisleyP]; %added one zero at end to make it match others

TScutoff = [1828,1828,1828]; %amount of time-series to cut off for spin up (5 years)

leapYearGuide=[365;365;365;366;365;365;365;366;365;365];
doy(:,1) = [1:1:leapYearGuide(1),1:1:leapYearGuide(2),1:1:leapYearGuide(3),...
        1:1:leapYearGuide(4),1:1:leapYearGuide(5),1:1:leapYearGuide(6),...
        1:1:leapYearGuide(7),1:1:leapYearGuide(8),1:1:leapYearGuide(9),1:1:leapYearGuide(10)];
doy(:,2) = [1:1:leapYearGuide(1),1:1:leapYearGuide(2),1:1:leapYearGuide(3),...
        1:1:leapYearGuide(4),1:1:leapYearGuide(5),1:1:leapYearGuide(6),...
        1:1:leapYearGuide(7),1:1:leapYearGuide(8),1:1:leapYearGuide(9),1:1:leapYearGuide(10)];  
leapYearGuide=[365;365;366;365;365;365;366;365;365;365];
doy(:,3) = [1:1:leapYearGuide(1),1:1:leapYearGuide(2),1:1:leapYearGuide(3),...
        1:1:leapYearGuide(4),1:1:leapYearGuide(5),1:1:leapYearGuide(6),...
        1:1:leapYearGuide(7),1:1:leapYearGuide(8),1:1:leapYearGuide(9),1:1:leapYearGuide(10)];

clearvars -except HySumclust options ndec doy rainTS TScutoff
        
%% Prescribed Water Fluxes
for EachRainTS = 1:3
P_t = round(rainTS(:,EachRainTS),ndec);
Pavg = mean(P_t); %average daily precipitation for Precip TS
ntri = length(HySumclust(:,1)); %number of montecarlo trials
%% Hydrologic Model Run
%QStsAll = cell(ntri,5); %pre-allocation
TTDbin = cell(ntri,11);
Qbin = cell(ntri,11);
for tri = 1:ntri
    %% KBM STORAGE - FLUX parameter choices
    eta = round(HySumclust(tri,1),ndec); %round to avoid truncation
    SuRef = round(HySumclust(tri,2),ndec);
    SlRef = round(HySumclust(tri,3),ndec);
    bu = round(HySumclust(tri,4),ndec);
    ku = (HySumclust(tri,5)); %dont round!
    bl = round(HySumclust(tri,6),ndec);
    kl = (HySumclust(tri,7)); %dont round!
    
    dt = 1; %timestep
    tlength = length(P_t); %simulation length (days)
    npart = 1*(10^(ndec)); %number of particles for age tracking
    %% KBM STORAGE - FLUX model initialization
    % Pre-allocation & initialization of storage & discharge vectors
    Sl_t = ones(tlength,1); %lower box storage
    Su_t = ones(tlength,1); %upper box storage
    Su_t(1) = SuRef; %initial set to reference storage upper box
    Sl_t(1) = SlRef; %initial set to reference storage lower box
    Ql_t = ones(tlength,2); %lower box discharge
    Qu_t = ones(tlength,3); %upper box discharge (total)
    Lr_t = ones(tlength,1); %upper box discharge (to lower box)
    % age tracking parameters 
    SuStorage = zeros(round(SuRef*npart),1);
    SlStorage = zeros(round(SlRef*npart),1);    
    %% KBM STORAGE - FLUX model Solution Scheme (with age track)
    ageDists = cell(length(P_t),1);
        for m = 2:length(P_t) % for each day of precipitation    
    % Previous time step values
            Sl = Sl_t(m-1); %lower storage
            Su = Su_t(m-1); %upper storage
            P =  P_t(m); % current precipitation (correct?)
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
    % Lower Box
	% stability check
            stcL = 0.5 + 0.5*(((Lr-kl*(Sl^bl))*dt)/( (Lr/kl)^(1/bl) - Sl));
            pl = min([stcL,1]);%stability criteria
            lowerBox = @(xl) (dt * ( Lr - pl*kl*(xl.^bl) - (1-pl)*kl*(Sl.^bl) )) - (xl - Sl);
            Sl_t(m) = fsolve(lowerBox,0,options); %lower box storage
            Ql_t(m,1) = Lr + ((Sl - Sl_t(m))/dt); %lower box outflow
            Ql_t(m,2) = Qu_t(m,1) +  Ql_t(m,1); %TOTAL discharge
    % Age track
            SuStorage = SuStorage+1; %update age in storage
            SlStorage = SlStorage+1; %update age in storage
            %add new water (preip) to upper box 
            SuStorage = [SuStorage;zeros(round(P*npart),1)]; 
            %randomly sample ages from upper box proportionally to outflow 
            remvUT = randsample(length(SuStorage),round(npart*L));
            %Save the exported ages
            outputAgesU = SuStorage(remvUT);
            %Remove exported ages from upper box
            SuStorage(remvUT) = NaN; SuStorage = SuStorage(isfinite(SuStorage));
            %Randomly sample exported ages to find which go to lower box
            remvUL = randsample(length(outputAgesU),round(npart*Lr));
            %Add outflow from upper box to lower box 
            SlStorage = [SlStorage;outputAgesU(remvUL)];
            %Remove ages added to lower box from upper box outflow
            %Remainder represents the age distribution added to discharge
            outputAgesU(remvUL) = NaN;  outputAgesU = outputAgesU(isfinite(outputAgesU));
            %randomly sample ages from lower box proportionally to outflow 
            remvLT = randsample(length(SlStorage),round(npart*Ql_t(m,1)));
            %Save the exported ages
            outputAgesL = SlStorage(remvLT);
            %Remove exported ages from lower box
            SlStorage(remvLT) = NaN; SlStorage = SlStorage(isfinite(SlStorage));
            %combine ages in discharge
            AgesInDischarge = [outputAgesU;outputAgesL];
            ageDists{m} = AgesInDischarge;
        end
    % Model Output Summary
    upperboxMB = [P_t,Qu_t(:,2),Su_t]; %UPPER input/total output/storage
    lowerboxMB = [Lr_t,Ql_t(:,1),Sl_t,Ql_t(:,2)]; %LOWER input/total output/storage/U+Ldischarge
    %% Extract Age Distributions
    startT = 1828; %start time to evaluate storage (i.e., cut initialization)
    Qw = lowerboxMB(startT:end,4); %total water discharge
    QBins = prctile(Qw,[10,20,30,40,50,60,70,80,90]); %percentiles to bin data
    aInd = find(Qw < QBins(1));
	bInd = find(Qw < QBins(2) & Qw >= QBins(1));
    cInd = find(Qw < QBins(3) & Qw >= QBins(2));
    dInd = find(Qw < QBins(4) & Qw >= QBins(3));
    eInd = find(Qw < QBins(5) & Qw >= QBins(4));
    fInd = find(Qw < QBins(6) & Qw >= QBins(5));
    gInd = find(Qw < QBins(7) & Qw >= QBins(6));
    hInd = find(Qw < QBins(8) & Qw >= QBins(7));
	iInd = find(Qw < QBins(9) & Qw >= QBins(8));
    jInd = find(Qw >= QBins(9));
    %Extract TTDs associated with each Q bin
    TTDall = ageDists(startT:end);
    TTDa = TTDall(aInd);
        TTDa = cat(1,TTDa{:});
	TTDb = TTDall(bInd);
        TTDb = cat(1,TTDb{:});
	TTDc = TTDall(cInd);
        TTDc = cat(1,TTDc{:});
	TTDd = TTDall(dInd);
        TTDd = cat(1,TTDd{:});
	TTDe = TTDall(eInd);
        TTDe = cat(1,TTDe{:});
	TTDf = TTDall(fInd);
        TTDf = cat(1,TTDf{:});
    TTDg = TTDall(gInd);
        TTDg = cat(1,TTDg{:});
	TTDh = TTDall(hInd);
        TTDh = cat(1,TTDh{:});
	TTDi = TTDall(iInd);
        TTDi = cat(1,TTDi{:});
	TTDj = TTDall(jInd);
        TTDj = cat(1,TTDj{:});
        
    TTDbin(tri,:) = {TTDall,TTDa,TTDb,TTDc,TTDd,TTDe,TTDf,TTDg,TTDh,TTDi,TTDj};
    Qbin(tri,:) = {Qw,Qw(aInd),Qw(bInd),Qw(cInd),Qw(dInd),Qw(eInd),...
        Qw(fInd),Qw(gInd),Qw(hInd),Qw(iInd),Qw(jInd)};     
end
save(strcat('QClustbin',num2str(EachRainTS),'.mat'),'Qbin')
save(strcat('TTDClustbin',num2str(EachRainTS),'.mat'),'TTDbin')
end