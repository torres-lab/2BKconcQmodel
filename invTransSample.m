function [dataSamp] = invTransSample(aData,nSamp)
[F,X] = ecdf(aData); %cdf of distribution
invCDF = fit(F,X,'linearinterp'); %inverse of CDF
dataSamp = feval(invCDF,rand(nSamp,1)); %random values from uniform dist


