function zdata = zscoreWholeMatrix(data)
% computes zscore using the mean & standard deviation of the whole dataset
dataSize = size(data); %find size of data matrix
dataList = reshape(data,[prod(dataSize),1]); %reshape to 1 column
dataListz = zscore(dataList); %zscore, 
zdata = reshape(dataListz,dataSize); %convert back to original shape






