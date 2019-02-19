function output = DSS_Analysis_Unpack(floc,filter,MLE)
%loads saved DSS analysis structure, unpacks structure, performs
%statistical computation and returns structure (excluding MLE estimates)

%% load data 

load(floc)
disp('Data loaded !')

%% Set whether post-processed data filtering is performed and whether maximum likelihoods are calculated
if nargin==1
    filter = 1;
    MLE = 0;
end

if nargin==2
    MLE = 0;
end

%% unpack input structure
disp('Unpacking data structures...')
DSS_Aux_structvars(output.input); 
DSS_Aux_structvars(output.Result);
%DSS_Aux_structvars(output.inputMC);
DSS_Aux_structvars(output.ResultMC);
disp('Data structures unpacked!')

%% include filter in ionisation rate

if filter
    disp('TeE data filter selected, filtering out erroneous TeE points, assuming a continuous decrease of TeE')
    
    if isfield(output.input, 'TeE_filterTime')
        tbeg = TeE_filterTime; 
        else
        tbeg=0.8; %beginning of TeE checks
    end
    %filter assumes at the end of the discharge a decrease in TeE
    TeE1MCS = TeE1MC;
    TeE2MCS = TeE2MC;
    
    for i=1:numel(n1Int(:,1)) %for each ROI
    %check if TeE1 is decreasing
        
    for j=1:numel(Iter)
    TeE1MCS(i,:,j) = DSS_Aux_smooth2a(squeeze(TeE1MC(i,:,j)),4); %smooth TeE to obtain rough TeE trend
    end
    
    V = diff(squeeze(TeE1MC(i,:,:)))>0; %Check where TeE increases
    V(output.input.Time < tbeg,:) = 0;
    V(2:end+1,:) = V;
    %turn results to NaN
    Ion1MC(i,V==1) = NaN; %Where it increases, push NaN into result
    PIon1MC(i,V==1) = NaN;
    TeE1MC(i,V==1) = NaN;
    CX1MC(i,V==1) = NaN;
    Rec1MC(i,V==1) = NaN;
    PRec1MC(i, V==1) = NaN;
    TeR1MC(i,V==1) = NaN;
    
    %check if TeE2 is decreasing
    for j=1:numel(Iter)
        TeE2MCS(i,:,j) = DSS_Aux_smooth2a(TeE2MC(i,:,j),4);
    end
    V = diff(squeeze(TeE2MC(i,:,:)))>0;
    V(output.input.Time < tbeg,:) = 0;
    V(2:end+1,:) = V;
    %use a safety buffer
    Ion2MC(i,V) = NaN; %Push NaN into result
    PIon2MC(i,V) = NaN;
    TeE2MC(i,V) = NaN;
    CX2MC(i,V) = NaN;
    Rec2MC(i,V) = NaN;
    PRec2MC(i,V) = NaN;
    TeR2MC(i,V) = NaN;
    
    end
    
end

disp('Filter done! Erroneous TeE entries overwritten with NaN')

%% Interpolate over missing data in Monte Carlo output

%initialise empty matrices for interpolated values across NaN
FRec1MCInan = 0.*FRec1MC+NaN;
FRec2MCInan = 0.*FRec2MC+NaN;
TeMCInan = 0.*TeMC+NaN;

IonMC1Inan = 0.*Ion1MC+NaN;
PIonMC1Inan = 0.*PIon1MC+NaN;
TeEMC1Inan = 0.*TeE1MC+NaN;
CXMC1Inan = 0.*CX1MC+NaN;
OK_Ion1 = sum(~isnan(Ion1MC),3)/Iter; %Entries with the fraction of non-filtered ionisation estimates from first Balmer line
output.Result.OK_Ion1 = OK_Ion1;

RecMC1Inan = 0.*Rec1MC+NaN;
TeRMC1Inan = 0.*TeR1MC+NaN;
PRecMC1Inan = 0.*PRec1MC+NaN;
OK_Rec1 = sum(~isnan(Rec1MC),3)/Iter; %Entries with the fraction of non-filtered recombination estimates from first Balmer line
output.Result.OK_Rec1 = OK_Rec1;

IonMC2Inan = 0.*Ion1MC+NaN;
PIonMC2Inan = 0.*PIon1MC+NaN;
TeEMC2Inan = 0.*TeE1MC+NaN;
CXMC2Inan = 0.*CX1MC+NaN;
OK_Ion2 = sum(~isnan(Ion2MC),3)/Iter; %Entries with the fraction of non-filtered ionisation estimates from second Balmer line
output.Result.OK_Ion2 = OK_Ion2;

RecMC2Inan = 0.*Rec1MC+NaN;
TeRMC2Inan = 0.*TeR1MC+NaN;
PRecMC2Inan = 0.*PRec1MC+NaN;
OK_Rec2 = sum(~isnan(Rec2MC),3)/Iter; %Entries with the fraction of non-filtered recombination estimates from second Balmer line
output.Result.OK_Rec2 = OK_Rec2;

n1Int(n1Int==0) = NaN;

disp('Interpolating NaN....')
%perform interpolation over NaN values (e.g. values which were filtered out
%or did not have a solution)
for i=1:numel(n1Int(:,1))
    for j=1:numel(n1Int(1,:))
        if ~isnan(n1Int(i,j))
            FRec1MCInan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(FRec1MC(i,j,:)),4);    
            FRec2MCInan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(FRec2MC(i,j,:)),4);
            TeMCInan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(TeMC(i,j,:)),4);

            IonMC1Inan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(Ion1MC(i,j,:)),4);
            PIonMC1Inan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(PIon1MC(i,j,:)),4);
            TeEMC1Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(TeE1MC(i,j,:)),4);
            CXMC1Inan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(CX1MC(i,j,:)), 4);
            RecMC1Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(Rec1MC(i,j,:)),4);
            TeRMC1Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(TeR1MC(i,j,:)),4);
            PRecMC1Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(PRec1MC(i,j,:)),4);

            IonMC2Inan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(Ion2MC(i,j,:)),4);
            PIonMC2Inan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(PIon2MC(i,j,:)),4);
            TeEMC2Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(TeE2MC(i,j,:)),4);
            CXMC2Inan(i,j,:) = DSS_Aux_inpaint_nans(squeeze(CX2MC(i,j,:)), 4);
            RecMC2Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(Rec2MC(i,j,:)),4);
            TeRMC2Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(TeR2MC(i,j,:)),4);
            PRecMC2Inan(i,j,: )= DSS_Aux_inpaint_nans(squeeze(PRec2MC(i,j,:)),4); 
        end    
    end
end
disp('Interpolation done!')

%% Perform integrals

disp('Performing integrals for integrated quantities...')
%initialise matrices
Rec1T = 0.*Time;
Ion1T = 0.*Time;
PRec1T = 0.*Time;
PIon1T = 0.*Time;
CX1T = 0.*Time;

%Initialise 
Rec1TMC = zeros(numel(Time),Iter) + NaN;
Ion1TMC = zeros(numel(Time),Iter) + NaN;
PRec1TMC = zeros(numel(Time),Iter) + NaN;
PIon1TMC = zeros(numel(Time),Iter) + NaN;
CX1TMC = zeros(numel(Time),Iter) + NaN;

%Initialise
Rec2T = 0.*Time;
Ion2T = 0.*Time;
PRec2T = 0.*Time;
PIon2T = 0.*Time;
CX2T = 0.*Time;

%Initialise
Rec2TMC = zeros(numel(Time),Iter) + NaN;
Ion2TMC = zeros(numel(Time),Iter) + NaN;
PRec2TMC = zeros(numel(Time),Iter) + NaN;
PIon2TMC = zeros(numel(Time),Iter) + NaN;
CX2TMC = zeros(numel(Time),Iter) + NaN;

%resolve NaN's in r,z position map
Z = DSS_Aux_inpaint_nans(Z,4);
R = DSS_Aux_inpaint_nans(R,4);

for i=1:numel(Time)
     
    if nansum(n1Int(:,i)./n1Int(:,i))>1 %If any non-NaN entries in Balmer line intensity
        
        %calculate integrals (single point - no MC) - first Balmer line 
        try %perform toroidal and poloidal integrals
        Rec1T(i) = Tor_Integr(Z(~isnan(Rec1(:,i)),i), R(~isnan(Rec1(:,i)),i), 2.*pi.*R(~isnan(Rec1(:,i)),i).*Rec1(~isnan(Rec1(:,i)),i)); 
        PIon1T(i) = Tor_Integr(Z(~isnan(PIon1(:,i)),i), R(~isnan(PIon1(:,i)),i), 2.*pi.*R(~isnan(PIon1(:,i)),i).*PIon1(~isnan(PIon1(:,i)),i));
        PRec1T(i) = Tor_Integr(Z(~isnan(PRec1(:,i)),i), R(~isnan(PRec1(:,i)),i), 2.*pi.*R(~isnan(PRec1(:,i)),i).*PRec1(~isnan(PRec1(:,i)),i));
        CX1T(i) = Tor_Integr(Z(~isnan(CX1(:,i)),i), R(~isnan(CX1(:,i)),i), 2.*pi.*R(~isnan(CX1(:,i)),i).*CX1(~isnan(CX1(:,i)),i));        
        Ion1T(i) = Tor_Integr(Z(~isnan(Ion1(:,i)),i), R(~isnan(Ion1(:,i)),i), 2.*pi.*R(~isnan(Ion1(:,i)),i).*Ion1(~isnan(Ion1(:,i)),i));
        catch
        end
    
        %calculate integrals (MC) - first Balmer line
        for j=1:Iter %perform toroidal and poloidal integrals
            Rec1TMC(i,j) = Tor_Integr(Z(~isnan(RecMC1Inan(:,i,j)),i), R(~isnan(RecMC1Inan(:,i,j)),i), 2.*pi.*R(~isnan(RecMC1Inan(:,i,j)),i).*RecMC1Inan(~isnan(RecMC1Inan(:,i,j)),i,j));
            Ion1TMC(i,j) = Tor_Integr(Z(~isnan(IonMC1Inan(:,i,j)),i),R(~isnan(IonMC1Inan(:,i,j)),i), 2.*pi.*R(~isnan(IonMC1Inan(:,i,j)),i).*IonMC1Inan(~isnan(IonMC1Inan(:,i,j)),i,j));
            PRec1TMC(i,j) = Tor_Integr(Z(~isnan(PRecMC1Inan(:,i,j)),i),R(~isnan(PRecMC1Inan(:,i,j)),i), 2.*pi.*R(~isnan(PRecMC1Inan(:,i,j)),i).*PRecMC1Inan(~isnan(PRecMC1Inan(:,i,j)),i,j));
            PIon1TMC(i,j) = Tor_Integr(Z(~isnan(PIonMC1Inan(:,i,j)),i),R(~isnan(PIonMC1Inan(:,i,j)),i), 2.*pi.*R(~isnan(PIonMC1Inan(:,i,j)),i).*PIonMC1Inan(~isnan(PIonMC1Inan(:,i,j)),i,j));
            CX1TMC(i,j) = Tor_Integr(Z(~isnan(CXMC1Inan(:,i,j)),i),R(~isnan(PIonMC1Inan(:,i,j)),i), 2.*pi.*R(~isnan(CXMC1Inan(:,i,j)),i).*CXMC1Inan(~isnan(CXMC1Inan(:,i,j)),i,j));        
        end
        
    end
    
    if nansum(n1Int(:,i)./n1Int(:,i))>1 %If any non-NaN entries in Balmer line intensity
        
        %calculate integrals (single point - no MC) - first Balmer line 
        try %perform toroidal and poloidal integrals
        Rec2T(i) = Tor_Integr(Z(~isnan(Rec2(:,i)),i), R(~isnan(Rec2(:,i)),i), 2.*pi.*R(~isnan(Rec2(:,i)),i).*Rec2(~isnan(Rec2(:,i)),i)); 
        PIon2T(i) = Tor_Integr(Z(~isnan(PIon2(:,i)),i), R(~isnan(PIon2(:,i)),i), 2.*pi.*R(~isnan(PIon2(:,i)),i).*PIon2(~isnan(PIon2(:,i)),i));
        PRec2T(i) = Tor_Integr(Z(~isnan(PRec2(:,i)),i), R(~isnan(PRec2(:,i)),i), 2.*pi.*R(~isnan(PRec2(:,i)),i).*PRec2(~isnan(PRec2(:,i)),i));
        CX2T(i) = Tor_Integr(Z(~isnan(CX2(:,i)),i), R(~isnan(CX2(:,i)),i), 2.*pi.*R(~isnan(CX2(:,i)),i).*CX2(~isnan(CX2(:,i)),i));        
        Ion2T(i) = Tor_Integr(Z(~isnan(Ion2(:,i)),i),R(~isnan(Ion2(:,i)),i),  2.*pi.*R(~isnan(Ion2(:,i)),i).*Ion2(~isnan(Ion2(:,i)),i));
        catch
        end
    
        %calculate integrals (MC) - first Balmer line
        for j=1:Iter %perform toroidal and poloidal integrals
            Rec2TMC(i,j) = Tor_Integr(Z(~isnan(RecMC2Inan(:,i,j)),i),R(~isnan(RecMC2Inan(:,i,j)),i), 2.*pi.*R(~isnan(RecMC2Inan(:,i,j)),i).*RecMC2Inan(~isnan(RecMC2Inan(:,i,j)),i,j));
            Ion2TMC(i,j) = Tor_Integr(Z(~isnan(IonMC2Inan(:,i,j)),i),R(~isnan(IonMC2Inan(:,i,j)),i), 2.*pi.*R(~isnan(IonMC2Inan(:,i,j)),i).*IonMC2Inan(~isnan(IonMC2Inan(:,i,j)),i,j));
            PRec2TMC(i,j) = Tor_Integr(Z(~isnan(PRecMC2Inan(:,i,j)),i),R(~isnan(PRecMC2Inan(:,i,j)),i), 2.*pi.*R(~isnan(PRecMC2Inan(:,i,j)),i).*PRecMC2Inan(~isnan(PRecMC2Inan(:,i,j)),i,j));
            PIon2TMC(i,j) = Tor_Integr(Z(~isnan(PIonMC2Inan(:,i,j)),i),R(~isnan(PIonMC2Inan(:,i,j)),i), 2.*pi.*R(~isnan(PIonMC2Inan(:,i,j)),i).*PIonMC2Inan(~isnan(PIonMC2Inan(:,i,j)),i,j));
            CX2TMC(i,j) = Tor_Integr(Z(~isnan(CXMC2Inan(:,i,j)),i),R(~isnan(CXMC2Inan(:,i,j)),i),2.*pi.*R(~isnan(CXMC2Inan(:,i,j)),i).*CXMC2Inan(~isnan(CXMC2Inan(:,i,j)),i,j));        
        end
        
    end   
    
end
disp('Integrations complete!')

%% Extract quantiles from output

disp('Extracting quantiles....')
FRec1Q = quantile(FRec1MCInan, Quant, 3); %Extracted quantiles from Monte Carlo distributions
FRec2Q = quantile(FRec2MCInan, Quant, 3);
TeQ = quantile(TeMCInan, Quant,3);
Ion1Q = quantile(IonMC1Inan, Quant, 3);
PIon1Q = quantile(PIonMC1Inan, Quant, 3);
TeE1Q = quantile(TeEMC1Inan, Quant, 3);
CX1Q = quantile(CXMC1Inan, Quant, 3);
CXIon1Q = quantile(CXMC1Inan./IonMC1Inan, Quant, 3);
Rec1Q = quantile(RecMC1Inan, Quant, 3);
TeR1Q = quantile(TeRMC1Inan, Quant, 3);
PRec1Q = quantile(PRecMC1Inan, Quant, 3);

Rec1TQ = quantile(Rec1TMC, Quant, 2); %Extracted quantiles from Monte Carlo distributions
Ion1TQ = quantile(Ion1TMC, Quant, 2);
PRec1TQ = quantile(PRec1TMC, Quant, 2);
PIon1TQ = quantile(PIon1TMC, Quant, 2);
CX1TQ = quantile(CX1TMC, Quant, 2);
CXIon1TQ = quantile(CX1TMC./Ion1TMC, Quant, 2);

Ion2Q = quantile(IonMC2Inan, Quant, 3); %Extracted quantiles from Monte Carlo distributions
PIon2Q = quantile(PIonMC2Inan, Quant, 3);
TeE2Q = quantile(TeEMC2Inan, Quant, 3);
CX2Q = quantile(CXMC2Inan, Quant, 3);
CXIon2Q = quantile(CXMC2Inan./IonMC2Inan, Quant, 3);
Rec2Q = quantile(RecMC2Inan, Quant, 3);
TeR2Q = quantile(TeRMC2Inan, Quant, 3);
PRec2Q = quantile(PRecMC2Inan, Quant, 3);

Rec2TQ = quantile(Rec2TMC, Quant, 2); %Extracted quantiles from Monte Carlo distributions
Ion2TQ = quantile(Ion2TMC, Quant, 2);
PRec2TQ = quantile(PRec2TMC, Quant, 2);
PIon2TQ = quantile(PIon2TMC, Quant, 2);
CX2TQ = quantile(CX2TMC, Quant, 2);
CXIon2TQ = quantile(CX2TMC./Ion2TMC, Quant, 2);
disp('Quantiles extracted!')

%% pack output variables in output structure

output.Result.FRec1Q = FRec1Q;
output.Result.FRec2Q = FRec2Q;
output.Result.TeQ = TeQ;
output.Result.Ion1Q = Ion1Q;
output.Result.Ion2Q = Ion2Q;
output.Result.PIon1Q = PIon1Q;
output.Result.PIon2Q = PIon2Q;
output.Result.TeE1Q = TeE1Q;
output.Result.TeE2Q = TeE2Q;
output.Result.CX1Q = CX1Q;
output.Result.CX2Q = CX2Q;
output.Result.CXIon1Q = CXIon1Q;
output.Result.CXIon2Q = CXIon2Q;
output.Result.Rec1Q = Rec1Q;
output.Result.Rec2Q = Rec2Q;
output.Result.PRec1Q = PRec1Q;
output.Result.PRec2Q = PRec2Q;
output.Result.TeR1Q = TeR1Q;
output.Result.TeR2Q = TeR2Q;

%Pack interpolated data
output.ResultMC.FRec1MCInan = FRec1MCInan;
output.ResultMC.FRec2MCInan = FRec2MCInan;
output.ResultMC.TeMCInan = TeMCInan;
output.ResultMC.IonMC1Inan = IonMC1Inan;
output.ResultMC.Ion1TMC = Ion1TMC;
output.ResultMC.IonMC2Inan = IonMC2Inan;
output.ResultMC.Ion2TMC = Ion2TMC;
output.ResultMC.PIonMC1Inan = PIonMC1Inan;
output.ResultMC.PIon1TMC = PIon1TMC;
output.ResultMC.PIonMC2Inan = PIonMC2Inan;
output.ResultMC.PIon2TMC = PIon2TMC;
output.ResultMC.TeEMC1Inan = TeEMC1Inan;
output.ResultMC.TeEMC2Inan = TeEMC2Inan;
output.ResultMC.CXMC1Inan = CXMC1Inan;
output.ResultMC.CX1TMC = CX1TMC;
output.ResultMC.CXMC2Inan = CXMC2Inan;
output.ResultMC.CX2TMC = CX2TMC;
output.ResultMC.RecMC1Inan = RecMC1Inan;
output.ResultMC.Rec1TMC = Rec1TMC;
output.ResultMC.RecMC2Inan = RecMC2Inan;
output.ResultMC.Rec2TMC = Rec2TMC;
output.ResultMC.PRecMC1Inan = PRecMC1Inan;
output.ResultMC.PRec1TMC = PRec1TMC;
output.ResultMC.PRecMC2Inan = PRecMC2Inan;
output.ResultMC.PRec2TMC = PRec2TMC;
output.ResultMC.TeRMC1Inan = TeRMC1Inan;
output.ResultMC.TeRMC2Inan = TeRMC2Inan;

output.Result.Ion1TQ = Ion1TQ;
output.Result.Ion2TQ = Ion2TQ;
output.Result.PIon1TQ = PIon1TQ;
output.Result.PIon2TQ = PIon2TQ;
output.Result.CX1TQ = CX1TQ;
output.Result.CX2TQ = CX2TQ;
output.Result.CXIon1TQ = CXIon1TQ;
output.Result.CXIon2TQ = CXIon2TQ;
output.Result.Rec1TQ = Rec1TQ;
output.Result.Rec2TQ = Rec2TQ;
output.Result.PRec1TQ = PRec1TQ;
output.Result.PRec2TQ = PRec2TQ;


%% Extract maximum likelihood from output

if MLE
    
Prob = 0.68; %68% confidence intervals

disp('Extracting Maximum Likelihood Estimates (MLE)....')
%initialise empty matrices
FRec1M = zeros(size(n1Int)) + NaN;
FRec2M = zeros(size(n1Int)) + NaN;
TeM = zeros(size(n1Int)) + NaN;
Ion1M = zeros(size(n1Int)) + NaN;
Ion2M = zeros(size(n1Int)) + NaN;
PIon1M = zeros(size(n1Int)) + NaN;
PIon2M = zeros(size(n1Int)) + NaN;
TeE1M = zeros(size(n1Int)) + NaN;
TeE2M = zeros(size(n1Int)) + NaN;
Rec1M = zeros(size(n1Int)) + NaN;
Rec2M = zeros(size(n1Int)) + NaN;
PRec1M = zeros(size(n1Int)) + NaN;
PRec2M = zeros(size(n1Int)) + NaN;
TeR1M = zeros(size(n1Int)) + NaN;
TeR2M = zeros(size(n1Int)) + NaN;
CX1M = zeros(size(n1Int)) + NaN;
CX2M = zeros(size(n1Int)) + NaN;
CXIon1M = zeros(size(n1Int)) + NaN;
CXIon2M = zeros(size(n1Int)) + NaN;
Ion1TM = zeros(size(Time)) + NaN;
Ion2TM = zeros(size(Time)) + NaN;
Rec1TM = zeros(size(Time)) + NaN;
Rec2TM = zeros(size(Time)) + NaN;
PIon1TM = zeros(size(Time)) + NaN;
PIon2TM = zeros(size(Time)) + NaN;
PRec1TM = zeros(size(Time)) + NaN;
PRec2TM = zeros(size(Time)) + NaN;
CX1TM = zeros(size(Time)) + NaN;
CX2TM = zeros(size(Time)) + NaN;
CXIon1TM = zeros(size(Time)) + NaN;
CXIon2TM = zeros(size(Time)) + NaN;

%Initialise lower limits
FRec1L = zeros(size(n1Int)) + NaN;
FRec2L = zeros(size(n1Int)) + NaN;
TeL = zeros(size(n1Int)) + NaN;
Ion1L = zeros(size(n1Int)) + NaN;
Ion2L = zeros(size(n1Int)) + NaN;
PIon1L = zeros(size(n1Int)) + NaN;
PIon2L = zeros(size(n1Int)) + NaN;
TeE1L = zeros(size(n1Int)) + NaN;
TeE2L = zeros(size(n1Int)) + NaN;
Rec1L = zeros(size(n1Int)) + NaN;
Rec2L = zeros(size(n1Int)) + NaN;
PRec1L = zeros(size(n1Int)) + NaN;
PRec2L = zeros(size(n1Int)) + NaN;
TeR1L = zeros(size(n1Int)) + NaN;
TeR2L = zeros(size(n1Int)) + NaN;
CX1L = zeros(size(n1Int)) + NaN;
CX2L = zeros(size(n1Int)) + NaN;
CXIon1L = zeros(size(n1Int)) + NaN;
CXIon2L = zeros(size(n1Int)) + NaN;
Ion1TL = zeros(size(Time)) + NaN;
Ion2TL = zeros(size(Time)) + NaN;
Rec1TL = zeros(size(Time)) + NaN;
Rec2TL = zeros(size(Time)) + NaN;
PIon1TL = zeros(size(Time)) + NaN;
PIon2TL = zeros(size(Time)) + NaN;
PRec1TL = zeros(size(Time)) + NaN;
PRec2TL = zeros(size(Time)) + NaN;
CX1TL = zeros(size(Time)) + NaN;
CX2TL = zeros(size(Time)) + NaN;
CXIon1TL = zeros(size(Time)) + NaN;
CXIon2TL = zeros(size(Time)) + NaN;

%Initialise upper limits
FRec1H = zeros(size(n1Int)) + NaN;
FRec2H = zeros(size(n1Int)) + NaN;
TeH = zeros(size(n1Int)) + NaN;
Ion1H = zeros(size(n1Int)) + NaN;
Ion2H = zeros(size(n1Int)) + NaN;
PIon1H = zeros(size(n1Int)) + NaN;
PIon2H = zeros(size(n1Int)) + NaN;
TeE1H = zeros(size(n1Int)) + NaN;
TeE2H = zeros(size(n1Int)) + NaN;
Rec1H = zeros(size(n1Int)) + NaN;
Rec2H = zeros(size(n1Int)) + NaN;
PRec1H = zeros(size(n1Int)) + NaN;
PRec2H = zeros(size(n1Int)) + NaN;
TeR1H = zeros(size(n1Int)) + NaN;
TeR2H = zeros(size(n1Int)) + NaN;
CX1H = zeros(size(n1Int)) + NaN;
CX2H = zeros(size(n1Int)) + NaN;
CXIon1H = zeros(size(n1Int)) + NaN;
CXIon2H = zeros(size(n1Int)) + NaN;
Ion1TH = zeros(size(Time)) + NaN;
Ion2TH = zeros(size(Time)) + NaN;
Rec1TH = zeros(size(Time)) + NaN;
Rec2TH = zeros(size(Time)) + NaN;
PIon1TH = zeros(size(Time)) + NaN;
PIon2TH = zeros(size(Time)) + NaN;
PRec1TH = zeros(size(Time)) + NaN;
PRec2TH = zeros(size(Time)) + NaN;
CX1TH = zeros(size(Time)) + NaN;
CX2TH = zeros(size(Time)) + NaN;
CXIon1TH = zeros(size(Time)) + NaN;
CXIon2TH = zeros(size(Time)) + NaN;

%smooth2a on Inan for less noisy MLEs

for i=1:numel(n1Int(:,1))
FRec1MCInan(i,:,:) = DSS_Aux_smooth2a(squeeze(FRec1MCInan(i,:,:)),0,1);
FRec2MCInan(i,:,:) = DSS_Aux_smooth2a(squeeze(FRec2MCInan(i,:,:)),0,1);
TeMCInan(i,:,:) = DSS_Aux_smooth2a(squeeze(TeMCInan(i,:,:)),0,1);
IonMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(IonMC1Inan(i,:,:)),0,1);
IonMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(IonMC2Inan(i,:,:)),0,1);
PIonMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(PIonMC1Inan(i,:,:)),0,1);
PIonMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(PIonMC2Inan(i,:,:)),0,1);
TeEMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(TeEMC1Inan(i,:,:)),0,1);
TeEMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(TeEMC2Inan(i,:,:)),0,1);
CXMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(CXMC1Inan(i,:,:)),0,1);
CXMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(CXMC2Inan(i,:,:)),0,1);
RecMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(RecMC1Inan(i,:,:)),0,1);
RecMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(RecMC2Inan(i,:,:)),0,1);
PRecMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(PRecMC1Inan(i,:,:)),0,1);
PRecMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(PRecMC2Inan(i,:,:)),0,1);
TeRMC1Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(TeRMC1Inan(i,:,:)),0,1);
TeRMC2Inan(i,:,:) = DSS_Aux_smooth2a(squeeze(TeRMC2Inan(i,:,:)),0,1);
end

%Again, smoothing
Ion1TMC = DSS_Aux_smooth2a(Ion1TMC,0,1);
Ion2TMC = DSS_Aux_smooth2a(Ion2TMC,0,1);
PIon1TMC = DSS_Aux_smooth2a(PIon1TMC,0,1);
PIon2TMC = DSS_Aux_smooth2a(PIon2TMC,0,1);
CX1TMC = DSS_Aux_smooth2a(CX1TMC,0,1);
CX2TMC = DSS_Aux_smooth2a(CX2TMC,0,1);
Rec1TMC = DSS_Aux_smooth2a(Rec1TMC,0,1);
Rec2TMC = DSS_Aux_smooth2a(Rec2TMC,0,1);
PRec1TMC = DSS_Aux_smooth2a(PRec1TMC,0,1);
PRec2TMC = DSS_Aux_smooth2a(PRec2TMC,0,1);

tic
%start MLE estimator calculation
for j=1:numel(n1Int(1,:))
    for i=1:numel(n1Int(:,1))
        if ~isnan(n1Int(i,j))
            i
            %Maximum likelihood estimation and uncertainty estimator
            %using DSS_MLE_Estimator
            [FRec1M(i,j), FRec1L(i,j), FRec1H(i,j)] = DSS_MLE_Estimator(squeeze(FRec1MCInan(i,j,:)), Prob);
            [FRec2M(i,j), FRec2L(i,j), FRec2H(i,j)] = DSS_MLE_Estimator(squeeze(FRec2MCInan(i,j,:)), Prob);
            [TeM(i,j), TeL(i,j), TeH(i,j)] = DSS_MLE_Estimator(squeeze(TeMCInan(i,j,:)), Prob);
            [Ion1M(i,j), Ion1L(i,j), Ion1H(i,j)] = DSS_MLE_Estimator(squeeze(IonMC1Inan(i,j,:)), Prob);
            [Ion2M(i,j), Ion2L(i,j), Ion2H(i,j)] = DSS_MLE_Estimator(squeeze(IonMC2Inan(i,j,:)), Prob);            
            [PIon1M(i,j), PIon1L(i,j), PIon1H(i,j)] = DSS_MLE_Estimator(squeeze(PIonMC1Inan(i,j,:)), Prob);
            [PIon2M(i,j), PIon2L(i,j), PIon2H(i,j)] = DSS_MLE_Estimator(squeeze(PIonMC2Inan(i,j,:)), Prob);                       
            [TeE1M(i,j), TeE1L(i,j), TeE1H(i,j)] = DSS_MLE_Estimator(squeeze(TeEMC1Inan(i,j,:)), Prob);
            [TeE2M(i,j), TeE2L(i,j), TeE2H(i,j)] = DSS_MLE_Estimator(squeeze(TeEMC2Inan(i,j,:)), Prob);                       
            [CX1M(i,j), CX1L(i,j), CX1H(i,j)] = DSS_MLE_Estimator(squeeze(CXMC1Inan(i,j,:)), Prob);
            [CX2M(i,j), CX2L(i,j), CX2H(i,j)] = DSS_MLE_Estimator(squeeze(CXMC2Inan(i,j,:)), Prob);                       
            [CXIon1M(i,j), CXIon1L(i,j), CXIon1H(i,j)] = DSS_MLE_Estimator(squeeze(CXMC1Inan(i,j,:)./IonMC1Inan(i,j,:)), Prob);
            [CXIon2M(i,j), CXIon2L(i,j), CXIon2H(i,j)] = DSS_MLE_Estimator(squeeze(CXMC2Inan(i,j,:)./IonMC2Inan(i,j,:)), Prob);                       
            [Rec1M(i,j), Rec1L(i,j), Rec1H(i,j)] = DSS_MLE_Estimator(squeeze(RecMC1Inan(i,j,:)), Prob);
            [Rec2M(i,j), Rec2L(i,j), Rec2H(i,j)] = DSS_MLE_Estimator(squeeze(RecMC2Inan(i,j,:)), Prob);            
            [PRec1M(i,j), PRec1L(i,j), PRec1H(i,j)] = DSS_MLE_Estimator(squeeze(PRecMC1Inan(i,j,:)), Prob);
            [PRec2M(i,j), PRec2L(i,j), PRec2H(i,j)] = DSS_MLE_Estimator(squeeze(PRecMC2Inan(i,j,:)), Prob);                       
            [TeR1M(i,j), TeR1L(i,j), TeR1H(i,j)] = DSS_MLE_Estimator(squeeze(TeRMC1Inan(i,j,:)), Prob);
            [TeR2M(i,j), TeR2L(i,j), TeR2H(i,j)] = DSS_MLE_Estimator(squeeze(TeRMC2Inan(i,j,:)), Prob);                       

%             [density,xmesh]=DSS_Aux_akde1d(squeeze(FRec1MCInan(i,j,:)));
%             [~,indx] = max(density);
%             FRec1M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(FRec2MCInan(i,j,:)));
%             [~,indx] = max(density);
%             FRec2M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(TeMCInan(i,j,:)));
%             [~,indx] = max(density);
%             TeM(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(IonMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             Ion1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(IonMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             Ion2M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PIonMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             PIon1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PIonMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             PIon2M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(TeEMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             TeE1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(TeEMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             TeE2M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CXMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             CX1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CXMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             CX2M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CXMC1Inan(i,j,:)./IonMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             CXIon1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CXMC2Inan(i,j,:)./IonMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             CXIon2M(i,j) = xmesh(indx);
%      
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(RecMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             Rec1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(RecMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             Rec2M(i,j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PRecMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             PRec1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PRecMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             PRec2M(i,j) = xmesh(indx);
%                        
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(TeRMC1Inan(i,j,:)));
%             [~,indx] = max(density);
%             TeR1M(i,j) = xmesh(indx);
%               
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(TeRMC2Inan(i,j,:)));
%             [~,indx] = max(density);
%             TeR2M(i,j) = xmesh(indx);

        end
    end
    %ML estimates non-profile parameter (e.g. integrated) parameters
        if nansum(n1Int(:,j)./n1Int(:,j))>1
            
            [Ion1TM(j), Ion1TL(j), Ion1TH(j)] = DSS_MLE_Estimator(squeeze(Ion1TMC(j,:)), Prob);
            [Ion2TM(j), Ion2TL(j), Ion2TH(j)] = DSS_MLE_Estimator(squeeze(Ion2TMC(j,:)), Prob);            
            [PIon1TM(j), PIon1TL(j), PIon1TH(j)] = DSS_MLE_Estimator(squeeze(PIon1TMC(j,:)), Prob);
            [PIon2TM(j), PIon2TL(j), PIon2TH(j)] = DSS_MLE_Estimator(squeeze(PIon2TMC(j,:)), Prob);                       
            [CX1TM(j), CX1TL(j), CX1TH(j)] = DSS_MLE_Estimator(squeeze(CX1TMC(j,:)), Prob);
            [CX2TM(j), CX2TL(j), CX2TH(j)] = DSS_MLE_Estimator(squeeze(CX2TMC(j,:)), Prob);                       
            [CXIon1TM(j), CXIon1TL(j), CXIon1TH(j)] = DSS_MLE_Estimator(squeeze(CX1TMC(j,:)./Ion1TMC(j,:)), Prob);
            [CXIon2TM(j), CXIon2TL(j), CXIon2TH(j)] = DSS_MLE_Estimator(squeeze(CX2TMC(j,:)./Ion2TMC(j,:)), Prob);                       
            [Rec1TM(j), Rec1TL(j), Rec1TH(j)] = DSS_MLE_Estimator(squeeze(Rec1TMC(j,:)), Prob);
            [Rec2TM(j), Rec2TL(j), Rec2TH(j)] = DSS_MLE_Estimator(squeeze(Rec2TMC(j,:)), Prob);            
            [PRec1TM(j), PRec1TL(j), PRec1TH(j)] = DSS_MLE_Estimator(squeeze(PRec1TMC(j,:)), Prob);
            [PRec2TM(j), PRec2TL(j), PRec2TH(j)] = DSS_MLE_Estimator(squeeze(PRec2TMC(j,:)), Prob);                       
             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(Ion1TMC(j,:)));
%             [~,indx] = max(density);
%             Ion1TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(Ion2TMC(j,:)));
%             [~,indx] = max(density);
%             Ion2TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(Rec1TMC(j,:)));
%             [~,indx] = max(density);
%             Rec1TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(Rec2TMC(j,:)));
%             [~,indx] = max(density);
%             Rec2TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PIon1TMC(j,:)));
%             [~,indx] = max(density);
%             PIon1TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PIon2TMC(j,:)));
%             [~,indx] = max(density);
%             PIon2TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PRec1TMC(j,:)));
%             [~,indx] = max(density);
%             PRec1TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(PRec2TMC(j,:)));
%             [~,indx] = max(density);
%             PRec2TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CX1TMC(j,:)));
%             [~,indx] = max(density);
%             CX1TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CX2TMC(j,:)));
%             [~,indx] = max(density);
%             CX2TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CX1TMC(j,:)./Ion1TMC(j,:)));
%             [~,indx] = max(density);
%             CXIon1TM(j) = xmesh(indx);
%             
%             [density,xmesh]=DSS_Aux_akde1d(squeeze(CX2TMC(j,:)./Ion2TMC(j,:)));
%             [~,indx] = max(density);
%             CXIon2TM(j) = xmesh(indx);
        end
        
        toc
        
end
disp('MLEs extracted!')

%% pack output MLE in output structure and save
output.Result.FRec1M = FRec1M;
output.Result.FRec2M = FRec2M;
output.Result.TeM = TeM;
output.Result.Ion1M = Ion1M;
output.Result.Ion2M = Ion2M;
output.Result.PIon1M = PIon1M;
output.Result.PIon2M = PIon2M;
output.Result.TeE1M = TeE1M;
output.Result.TeE2M = TeE2M;
output.Result.CX1M = CX1M;
output.Result.CX2M = CX2M;
output.Result.CXIon1M = CXIon1M;
output.Result.CXIon2M = CXIon2M;
output.Result.Rec1M = Rec1M;
output.Result.Rec2M = Rec2M;
output.Result.PRec1M = PRec1M;
output.Result.PRec2M = PRec2M;
output.Result.TeR1M = TeR1M;
output.Result.TeR2M = TeR2M;

output.Result.Ion1TM = Ion1TM;
output.Result.Ion2TM = Ion2TM;
output.Result.PIon1TM = PIon1TM;
output.Result.PIon2TM = PIon2TM;
output.Result.CX1TM = CX1TM;
output.Result.CX2TM = CX2TM;
output.Result.CXIon1TM = CXIon1TM;
output.Result.CXIon2TM = CXIon2TM;
output.Result.Rec1TM = Rec1TM;
output.Result.Rec2TM = Rec2TM;
output.Result.PRec1TM = PRec1TM;
output.Result.PRec2TM = PRec2TM;

output.Result.FRec1L = FRec1L;
output.Result.FRec2L = FRec2L;
output.Result.TeL = TeL;
output.Result.Ion1L = Ion1L;
output.Result.Ion2L = Ion2L;
output.Result.PIon1L = PIon1L;
output.Result.PIon2L = PIon2L;
output.Result.TeE1L = TeE1L;
output.Result.TeE2L = TeE2L;
output.Result.CX1L = CX1L;
output.Result.CX2L = CX2L;
output.Result.CXIon1L = CXIon1L;
output.Result.CXIon2L = CXIon2L;
output.Result.Rec1L = Rec1L;
output.Result.Rec2L = Rec2L;
output.Result.PRec1L = PRec1L;
output.Result.PRec2L = PRec2L;
output.Result.TeR1L = TeR1L;
output.Result.TeR2L = TeR2L;

output.Result.Ion1TL = Ion1TL;
output.Result.Ion2TL = Ion2TL;
output.Result.PIon1TL = PIon1TL;
output.Result.PIon2TL = PIon2TL;
output.Result.CX1TL = CX1TL;
output.Result.CX2TL = CX2TL;
output.Result.CXIon1TL = CXIon1TL;
output.Result.CXIon2TL = CXIon2TL;
output.Result.Rec1TL = Rec1TL;
output.Result.Rec2TL = Rec2TL;
output.Result.PRec1TL = PRec1TL;
output.Result.PRec2TL = PRec2TL;

output.Result.FRec1H = FRec1H;
output.Result.FRec2H = FRec2H;
output.Result.TeH = TeH;
output.Result.Ion1H = Ion1H;
output.Result.Ion2H = Ion2H;
output.Result.PIon1H = PIon1H;
output.Result.PIon2H = PIon2H;
output.Result.TeE1H = TeE1H;
output.Result.TeE2H = TeE2H;
output.Result.CX1H = CX1H;
output.Result.CX2H = CX2H;
output.Result.CXIon1H = CXIon1H;
output.Result.CXIon2H = CXIon2H;
output.Result.Rec1H = Rec1H;
output.Result.Rec2H = Rec2H;
output.Result.PRec1H = PRec1H;
output.Result.PRec2H = PRec2H;
output.Result.TeR1H = TeR1H;
output.Result.TeR2H = TeR2H;

output.Result.Ion1TH = Ion1TH;
output.Result.Ion2TH = Ion2TH;
output.Result.PIon1TH = PIon1TH;
output.Result.PIon2TH = PIon2TH;
output.Result.CX1TH = CX1TH;
output.Result.CX2TH = CX2TH;
output.Result.CXIon1TH = CXIon1TH;
output.Result.CXIon2TH = CXIon2TH;
output.Result.Rec1TH = Rec1TH;
output.Result.Rec2TH = Rec2TH;
output.Result.PRec1TH = PRec1TH;
output.Result.PRec2TH = PRec2TH;

disp('Saving MLEs...') %Save MLEs

save(floc, 'output','-append')

%DSS_Analysis_FileCleanup(shotnr,linecombo)

end

if NOLAC == 1
    disp('File does not belong on LAC, removing')
    delete(floc)
end

disp('Done. Returning output structure')
end

function Z=Tor_Integr(X,Y,Integr)
%Function for toroidal/poloidal integral

[~,Xi] = arclength(X,Y,'spline'); 
Xi = [0; Xi;]; %Assemble a distance vector along the path given by the x,y positions
Z = trapz(cumsum(Xi),Integr); %Integrate the field Integr along this distance vector

end

function [arclen,seglen] = arclength(px,py,varargin)
% arclength: compute arc length of a space curve, or any curve represented as a sequence of points
% usage: [arclen,seglen] = arclength(px,py)         % a 2-d curve
% usage: [arclen,seglen] = arclength(px,py,pz)      % a 3-d space curve
% usage: [arclen,seglen] = arclength(px,py,method)  % specifies the method used
%
% Computes the arc length of a function or any
% general 2-d, 3-d or higher dimensional space
% curve using various methods.
%
% arguments: (input)
%  px, py, pz, ... - vectors of length n, defining points
%        along the curve. n must be at least 2. Replicate
%        points should not be present in the curve.
%
%  method - (OPTIONAL) string flag - denotes the method
%        used to compute the arc length of the curve.
%
%        method may be any of 'linear', 'spline', or 'pchip',
%        or any simple contraction thereof, such as 'lin',
%        'sp', or even 'p'.
%        
%        method == 'linear' --> Uses a linear chordal
%               approximation to compute the arc length.
%               This method is the most efficient.
%
%        method == 'pchip' --> Fits a parametric pchip
%               approximation, then integrates the
%               segments numerically.
%
%        method == 'spline' --> Uses a parametric spline
%               approximation to fit the curves, then
%               integrates the segments numerically.
%               Generally for a smooth curve, this
%               method may be most accurate.
%
%        DEFAULT: 'linear'
%
%
% arguments: (output)
%  arclen - scalar total arclength of all curve segments
%
%  seglen - arclength of each independent curve segment
%           there will be n-1 segments for which the
%           arc length will be computed.
%
%
% Example:
% % Compute the length of the perimeter of a unit circle
% theta = linspace(0,2*pi,10);
% x = cos(theta);
% y = sin(theta);
%
% % The exact value is
% 2*pi
% % ans =
% %          6.28318530717959
%
% % linear chord lengths
% arclen = arclength(x,y,'l')
% % arclen =
% %           6.1564
%
% % Integrated pchip curve fit
% arclen = arclength(x,y,'p')
% % arclen =
% %          6.2782
%
% % Integrated spline fit
% arclen = arclength(x,y,'s')
% % arclen =
% %           6.2856
%
% Example:
% % A (linear) space curve in 5 dimensions
% x = 0:.25:1;
% y = x;
% z = x;
% u = x;
% v = x;
%
% % The length of this curve is simply sqrt(5)
% % since the "curve" is merely the diagonal of a
% % unit 5 dimensional hyper-cube.
% [arclen,seglen] = arclength(x,y,z,u,v,'l')
% % arclen =
% %           2.23606797749979
% % seglen =
% %         0.559016994374947
% %         0.559016994374947
% %         0.559016994374947
% %         0.559016994374947
%
%
% See also: interparc, spline, pchip, interp1
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release: 1.0
% Release date: 3/10/2010

% unpack the arguments and check for errors
if nargin < 2
  error('ARCLENGTH:insufficientarguments', ...
    'at least px and py must be supplied')
end

n = length(px);
% are px and py both vectors of the same length?
if ~isvector(px) || ~isvector(py) || (length(py) ~= n)
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of the same length')
elseif n < 2
  error('ARCLENGTH:improperpxorpy', ...
    'px and py must be vectors of length at least 2')
end

% compile the curve into one array
data = [px(:),py(:)];

% defaults for method and tol
method = 'linear';

% which other arguments are included in varargin?
if numel(varargin) > 0
  % at least one other argument was supplied
  for i = 1:numel(varargin)
    arg = varargin{i};
    if ischar(arg)
      % it must be the method
      validmethods = {'linear' 'pchip' 'spline'};
      ind = strmatch(lower(arg),validmethods);
      if isempty(ind) || (length(ind) > 1)
        error('ARCLENGTH:invalidmethod', ...
          'Invalid method indicated. Only ''linear'',''pchip'',''spline'' allowed.')
      end
      method = validmethods{ind};
      
    else
      % it must be pz, defining a space curve in higher dimensions
      if numel(arg) ~= n
        error('ARCLENGTH:inconsistentpz', ...
          'pz was supplied, but is inconsistent in size with px and py')
      end
      
      % expand the data array to be a 3-d space curve
      data = [data,arg(:)]; %#ok
    end
  end
  
end

% what dimension do we live in?
nd = size(data,2);

% compute the chordal linear arclengths
seglen = sqrt(sum(diff(data,[],1).^2,2));
arclen = sum(seglen);

% we can quit if the method was 'linear'.
if strcmpi(method,'linear')
  % we are now done. just exit
  return
end

% 'spline' or 'pchip' must have been indicated,
% so we will be doing an integration. Save the
% linear chord lengths for later use.
chordlen = seglen;

% compute the splines
spl = cell(1,nd);
spld = spl;
diffarray = [3 0 0;0 2 0;0 0 1;0 0 0];
for i = 1:nd
  switch method
    case 'pchip'
      spl{i} = pchip([0;cumsum(chordlen)],data(:,i));
    case 'spline'
      spl{i} = spline([0;cumsum(chordlen)],data(:,i));
      nc = numel(spl{i}.coefs);
      if nc < 4
        % just pretend it has cubic segments
        spl{i}.coefs = [zeros(1,4-nc),spl{i}.coefs];
        spl{i}.order = 4;
      end
  end
  
  % and now differentiate them
  xp = spl{i};
  xp.coefs = xp.coefs*diffarray;
  xp.order = 3;
  spld{i} = xp;
end

% numerical integration along the curve
polyarray = zeros(nd,3);
for i = 1:spl{1}.pieces
  % extract polynomials for the derivatives
  for j = 1:nd
    polyarray(j,:) = spld{j}.coefs(i,:);
  end
  
  % integrate the arclength for the i'th segment
  % using quadgk for the integral. I could have
  % done this part with an ode solver too.
  seglen(i) = quadgk(@(t) segkernel(t),0,chordlen(i));
end

% and sum the segments
arclen = sum(seglen);

% ==========================
%   end main function
% ==========================
%   begin nested functions
% ==========================
  function val = segkernel(t)
    % sqrt((dx/dt)^2 + (dy/dt)^2)
    
    val = zeros(size(t));
    for k = 1:nd
      val = val + polyval(polyarray(k,:),t).^2;
    end
    val = sqrt(val);
    
  end % function segkernel

end % function arclength



