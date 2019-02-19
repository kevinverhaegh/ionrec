function [MLE, Q1, Q2] = DSS_MLE_Estimator(Y, Prob) 
%Function which, given a distribution of points, maps a probability density
%function, locates the peak (e.g. Maximum Likelihood) and determines the
%smallest possible interval over which the integral of the PDF amounts to
%the set probability - e.g. Highest Density Interval

if nargin==1
    Prob = 0.68; %default 68% confidence interval
end

Y = Y(:); %re-organise data and remove NaN
Y = Y(~isnan(Y));

[density,xmesh]=DSS_Aux_akde1d(Y); %Determine PDF using Botev's algorithm
if ~isnan(density) %If PDF extraction is successfull

Val = cumtrapz(xmesh, density); %Cumulitive integral

Interv = [];
Q1 = [];
Q2 = [];


for i=1:(numel(xmesh)-1) %Loop over x-mesh
    [~,indx] = min(abs(Val - Val(i) - Prob)); %Find where the difference in integrated probability to the i-point is Prob
    if abs(Val(indx) - Val(i) - Prob) < 0.005 %If the difference between both points is less than a set tolerance
        Interv(end+1) = abs(xmesh(indx) - xmesh(i)); %Save the interval
        Q1(end+1) = xmesh(i); %Save the domain
        Q2(end+1) = xmesh(indx);
    end
end

[~,indx] = min(Interv); %Get the minimum interval
if isempty(indx)
    MLE = NaN;
    Q1 = NaN;
    Q2 = NaN;
else %Get Lower/Upper bounds from minimum interval
    Q1 = Q1(indx);
    Q2 = Q2(indx);
    xmeshAdj = xmesh(xmesh > Q1 & xmesh < Q2); %Adjust x-mesh axis
    densityAdj = density(xmesh > Q1 & xmesh < Q2); %Adjust PDF axis
    [~,indx] = max(densityAdj); %Get MLE position
    MLE = xmeshAdj(indx); %Store MLE
end
else
    MLE = NaN;
    Q1 = NaN;
    Q2 = NaN;
end