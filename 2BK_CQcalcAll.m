%% Parameter Choices
DaRange = [0.1 1 10 100];
LoadData = 0; % do you need to load saved TTD and Q data? 1 = yes, 0 = no
%%
if LoadData == 1
 %% Load Data
load('TTDClustbin1.mat');
SRttds = TTDbin;
 load('TTDClustbin2.mat');
BRttds = TTDbin;
 load('TTDClustbin3.mat');
PRttds = TTDbin;
% 
 load('QClustbin1.mat');
SRQws = Qbin;
 load('QClustbin2.mat');
BRQws = Qbin;
 load('QClustbin3.mat');
PRQws = Qbin;
%
clear TTDbin Qbin
% Combine TTD and Q Data into one matrix
TTDbin = [SRttds(:,2:end);BRttds(:,2:end);PRttds(:,2:end);];
Qbin = [SRQws(:,2:end);BRQws(:,2:end);PRQws(:,2:end);];
%% Calculate mean Q and transit T
Qmeans = zeros(90,10); %mean Q for each bin
TTDmean = zeros(90,10); %mean TTD for each bin
TTDmed = zeros(90,10); %median TTD for each bin
TTDtmean = zeros(90,10); %trimmed mean TTD for each bin
for m = 1:length(Qmeans(:,1))
    for mm = 1:length(Qmeans(1,:))
            Qmeans(m,mm) = mean(Qbin{m,mm});
			TTDfromBin = TTDbin{m,mm};
			TTDmean(m,mm) = mean(TTDfromBin);
			TTDmed(m,mm) = median(TTDfromBin);
			trimRange = prctile(	TTDfromBin,[25,75]);
			trimRangeIDX = find(TTDfromBin>=trimRange(1) & TTDfromBin<=trimRange(2));	
			TTDtmean(m,mm) = mean(TTDfromBin(trimRangeIDX));
    end 
end
wmeanTT = wmean(TTDmean,Qmeans,2); %mean TT weighted by mean Q
wmedTT = wmean(TTDmed,Qmeans,2); %median TT weighted by mean Q
wtmeanTT = wmean(TTDtmean,Qmeans,2); %trimmed mean TT weighted by mean Q
%Load Kinetic Model Output
load('kineticParam.mat');
load('clustModels.mat');
load('kineticFit.mat');
else
end

%% Full dataset C-Q analysis
bValsNa = zeros(90,20,4); %pre-allocation (90 X 20 X 4)
bValsSi = zeros(90,20,4); %pre-allocation
% Set Kinetic Parameters
for WModel = 1:20 %pick weathering model
    NaWfunc = @(t) 1 - ((1-kineticParam(clustModels(WModel),4)).*...
        (kineticParam(clustModels(WModel),8)-t)).^...
        (1./(1-kineticParam(clustModels(WModel),4)));
    SiWfunc = kineticFit{clustModels(WModel),2}; %interp of Si model
    WModelTmax = kineticParam(clustModels(WModel),5);
    for DaRangeIdx = 1:4 %set Da
        Da = DaRange(DaRangeIdx); %[0.1 1 10 100] set Da
        %% Select hydroparameters
        parfor hyparams = 1:90 %pick hydrologic model
            %% Open discharge bins
            cVals = zeros(10,2);
            for binQ = 1:10	% for each discharge bin
                td = (TTDbin{hyparams,binQ}./wmeanTT(hyparams)); %norm TTD
                tdRS = Da.* invTransSample(td,1E6); %samples from TTD
                concs = zeros(1E6,2); %conc pre-allocation
                concs(tdRS>=WModelTmax,2) = 1; %if TT is long, conc = 1
                %solve for conc
                concs(:,1) = NaWfunc(tdRS);
                concs(tdRS<WModelTmax,2) = feval(SiWfunc,tdRS(tdRS<WModelTmax));
                cVals(binQ,1) = mean(concs(:,1)); %take mean of conc. dist
                cVals(binQ,2) = mean(concs(:,2)); %take mean of conc. dist
            end
        %% Calculate C-Q power law fits
        pwrFitANa = fit(Qmeans(hyparams,:)',cVals(:,1),'power1');
        pwrFitAcoefNa = coeffvalues(pwrFitANa);
        pwrFitASi = fit(Qmeans(hyparams,:)',cVals(:,2),'power1');
        pwrFitAcoefSi = coeffvalues(pwrFitASi);
        %% Save Data
        bValsNa(hyparams,WModel,DaRangeIdx) = pwrFitAcoefNa(2);
        bValsSi(hyparams,WModel,DaRangeIdx) = pwrFitAcoefSi(2);
        end
    end
end

%% Reshape 3D matrices into single columns
dataD01(:,1) = reshape(bValsNa(:,:,1),[90*20,1]);
dataD01(:,2) = reshape(bValsSi(:,:,1),[90*20,1]);
dataD1(:,1) = reshape(bValsNa(:,:,2),[90*20,1]);
dataD1(:,2) = reshape(bValsSi(:,:,2),[90*20,1]);
dataD10(:,1) = reshape(bValsNa(:,:,3),[90*20,1]);
dataD10(:,2) = reshape(bValsSi(:,:,3),[90*20,1]);
dataD100(:,1) = reshape(bValsNa(:,:,4),[90*20,1]);
dataD100(:,2) = reshape(bValsSi(:,:,4),[90*20,1]);
% group all b values together into single structure
dataDall = [dataD01;dataD1;dataD10;dataD100];
%% Simple Plotting
DaGroup = [ones(1800,1).*0.1;ones(1800,1);ones(1800,1).*10;ones(1800,1).*100];
figure(3)
clf
scatterhist(dataDall(:,1),dataDall(:,2),'Group',DaGroup,'Kernel','on',...
    'direction','out')
xlabel('b-exponent Na'); ylabel('b-exponent Si')

