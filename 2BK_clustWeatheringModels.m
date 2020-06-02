%% Monte-Carlo Search
[kineticParam, kineticFit, dataVals] = WeatheringModelMCsample(1000);

%% Find Cluster end-members
nClust = 20;
[IDX,C] = kmeans(dataVals(:,:),nClust); %Cluster
clustModels = zeros(nClust,1); %20 clusters
for m = 1:nClust %for each cluseter
    clustInd = find(IDX == m); %select cluster
    diffCent = zeros(length(clustInd),1); %pre-allocation
    for mm = 1:length(clustInd) %for all models within cluster
        %difference model from the cluster centroid for Na and Si
        diffCent(mm,1) = sum((C(m,:) - dataVals(clustInd(mm),:)).^2);
    end
	zDif = find(diffCent(:,1) == min(diffCent(:,1))); %closest to centroid 
    clustModels(m,1) = clustInd(zDif); 
end
%% Save Data
save('kineticParam.mat',kineticParam);
save('clustModels.mat',clustModels);
save('kineticFit.mat',kineticFit);