function [kineticParam, kineticFit, dataVals] = WeatheringModelMCsample(ntrials)
warning('off','all')
%% Monte Carlo Parameter Search of Weathering Model
kineticParam = zeros(ntrials,8); %A B mp md Tend Teq Siend 
kineticFit = cell(ntrials,2);
dataVals = zeros(ntrials,200);
for tri = 1:ntrials
%% Random Parameter values
    mp = random('uniform',1,2); % precipitation rate law exponent
    md = random('uniform',1,2); % dissolution rate law exponent
    A = random('uniform',1,10); % NaEq./(SiEq.*u);
    B = 10.^(random('uniform',-1,2)); % Rp./Rd;
%% Fixed Parameters
    options = odeset('reltol',1e-5,'abstol',1e-9); 
    x0 = [0,0]; %initial values for diff. eq. solver
    t_end = 10; %minimum model run time
    SiEnd = 10; %initial condition for while loop
%% Diff EQ Solver
    %while loop to run model for as short as possible for Si to go to ~1
    while SiEnd > 1.01 %run model until Si reaches at least 1.01
        [t,x] = ode23tb(@weatheringModelEquations,[0,t_end],x0,options,A,B,md,mp);
        SiEnd = x(end,2); %Si at end of model
        if isreal(x(:,2)) == 0 %where there any imaginary solutions?
            t_end = t_end - 1; %if so, go back one timestep and re-run model
            [t,x] = ode23tb(@weatheringModelEquations,[0,t_end],x0,options,A,B,md,mp);
            break %exit loop
        else
            t_end = t_end+1; %run model for one timestep longer
        end
    end
    %analytical solution to Na conc.
    intConst = (((1-x(end,1)).^(1-md))./(1-md))+t(end); %integration constant
    Teq = -1.*((((1-0.99).^(1-md))./(1-md))-intConst); % t when Na* = 0.99
	interpFitSi = fit(t,x(:,2),'linearinterp'); %linear interpolation of Si
	interpFitNa = fit(t,x(:,1),'linearinterp'); %linear interpolation of Na
%% Save output 
    kineticFit{tri,1} = interpFitNa; %save fits in cell
	kineticFit{tri,2} = interpFitSi;
    xRange = linspace(0,11); %t range to evaluate fits
    dataVals(tri,1:100) = feval(interpFitNa,xRange); %evalaute fits for clustering
    dataVals(tri,101:200) = feval(interpFitSi,xRange);
    kineticParam(tri,:) = [A B mp md t_end Teq x(end,2) intConst]; %save parameters
end