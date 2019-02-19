function [SCD, ACD, PECrec, PECexc, PLT, PRB, CCD] = TECPEC(N, Ne, Te, ADS,indx)

%% Code overview

%This function determines the photon emission coefficients (m^3 / s) of recombination
%and excitation  of Balmer lines given a set of electron densities Ne (m^-3) and a
%set of electron temperatures Te (eV) given the deuterium ADAS data. Also
%returned are the neutral fP and ion fractional abundances fM (see assumptions
%made below). Note that pairs of Ne, Te need to be provided.

%June 2015 - Code creation - v 0.1 
%March 2016 - Code revision - v 0.2
%Oktober 2017 - Add inx parameter - v 0.3

%Kevin Verhaegh v 0.3 University of York / EPFL-SPC

%% Load deuterium adas and explanation on data set

if nargin>3
    AdasNe=ADS.AdasNe;
    AdasTe=ADS.AdasTe;
    SCDHydrogen=ADS.SCDHydrogen;
    ACDHydrogen=ADS.ACDHydrogen;
    CCDHydrogen = ADS.CCDHydrogen;
    BalmerRecombination=ADS.BalmerRecombination;
    BalmerExcitation=ADS.BalmerExcitation;
    PLTHydrogen = ADS.PLTHydrogen;
    PRBHydrogen = ADS.PRBHydrogen;
else
    
load('/home/verhaegh/DeuteriumADAS') %load ADAS data set

end

if nargin~=4
    indx = 1:7;
end

SCD=NaN;
ACD=NaN;
PECrec=NaN;
PECexc=NaN;
PLT=NaN;
PRB=NaN; 
CCD=NaN;

%DeuteriumADAS is a matlab dataset with data exported from:

% Open ADAS - Full Path - adf15/pec12#h/pec12#h_pju#h0.dat 
%(contains PECs for all transitions up until n=12 as function of ne, Te)

% Open ADAS - FUll Path - adf11/acd12/acd12_h.dat
%(contains effective recombination coefficients as function of ne, Te)
%in units cm^3 /s (note the original dataset is given in 10^{parameter
%noted}. This is not the case in the loaded file.

%Open ADAS - Full Path - adf11/scd12/scd12_h.dat
%(contains effective ionisation coefficients as function of ne, Te)
%in units cm^3 /s (note the original dataset is given in 10^{parameter
%noted}. This is not the case in the loaded file. 

%Note on ADAS effective ionisation/recombination coefficients:
%The recombination / ionisation rate provided by ADAS applies in an
%infinite plasma with homogeneous temperature and density. The balance of
%recombination / ionisation rate is often used to determine the fractional
%abundance. Note that this ignores all transport effects and provides the
%fractional abundance you would assume in an infinite plasma in thermal
%equilibrium with homogeneous temperature and density.

%DeuteriumADAS is ADAS data containing:
%ADAS Ne, Te, grid (AdasNe (m^{-3}), AdasTe (eV))
%ACDHydrogen; SCDHydrogen ((**respectively effective recombination and ionisation 
% coefficients in cm^{3} /s in [Te, ne])
%BalmerExcitation (PECs for excitation [transition, Te, ne]
%BalmerExcitationTn2 (upper n levels of Balmer transitions)
%BalmerRecombination (PECs for recombination [transition, Te, ne]
%lambdaBalmer.... (wavelength of Balmer transitions)

%% Initialise parameters
N = N(:);
Te = Te(:);
Ne = Ne(:);
PECrec = zeros(numel(N), numel(Ne));
PECexc = 0.*PECrec;

%% For all given densities

              %Fractional abundance calculation based on assuming an infinite
%             %plasma with homogeneous density/temperature and using CR
%             %effective recombination / ionisation rates (ignoring
%             %transport; assuming fractional abundance is only determined by
%             %a balance of effective ionisation / recombination rates)
if sum(find(indx==1))>0
SCD = 1e-6 * interp2L(AdasNe, AdasTe, SCDHydrogen', Ne, Te);
end
if sum(find(indx==2))>0
ACD = 1e-6 * interp2L(AdasNe, AdasTe, ACDHydrogen', Ne, Te);
end
if sum(find(indx==7))>0
CCD = 1e-6 * interp2L(AdasNe, AdasTe, CCDHydrogen', Ne, Te);
end
if sum(find(indx==6))>0
PRB = 1e-6 * interp2L(AdasNe, AdasTe, PRBHydrogen', Ne, Te);
%correct PRB for Bremsstrahlung
PRB = PRB - 5.35e-37.*sqrt(Te/1e3); %Wesson
PRB(PRB<0) = 0;
end
if sum(find(indx==5))>0
PLT = 1e-6 * interp2L(AdasNe, AdasTe, PLTHydrogen', Ne, Te);
end

if sum(find(indx==3))>0
PECrec = zeros(numel(N), numel(Ne));
for i=1:numel(N)
     PECrec(i,:) = 1e-6 * interp2L(AdasNe, AdasTe, reshape(BalmerRecombination(N(i),:,:), [29,24]), Ne, Te); %Converts the PEC from cm^3 / s to m^3 / s
end
end

if sum(find(indx==4))>0
PECexc = zeros(numel(N), numel(Ne));  
for i=1:numel(N)
PECexc(i,:) = 1e-6 * interp2L(AdasNe, AdasTe, reshape(BalmerExcitation(N(i),:,:), [29,24]), Ne, Te); %Converts the PEC from cm^3 / s to m^3 / s
end
end


function alpha = interp2L(x,y,V,xN,yN)

%performes 2D surface interpolation (spline) on a logaritmic grid

F = griddedInterpolant({log10(x), log10(y)}, log10(V'), 'linear','none');
alpha = 10.^F(log10(xN), log10(yN));
%alpha = 10.^interp2(log10(x), log10(y), log10(V), log10(xN), log10(yN),'linear');

