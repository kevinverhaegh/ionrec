function output = DSS_Analysis_Calc_RecIon(input)
%Function to determine Recombination / Ionisation rates

%% Analyse input, get necessary info from input structure

n1 = input.n1; %Balmer line index

IR = input.IR; %Balmer line intensity (recombinative emission
IE = input.IE; %Balmer line intensity (ionisation emission
none = input.none; %neutral fraction
Den = input.Den; %electron density
DL = input.DL; %Divertor leg length

%% Get optional info from input structure
if isfield(input, 'Te') %initial Te vector
    Te = input.Te;
else
    Te = 0.2:0.02:200; 
end

if isfield(input, 'ADASInfo') %ADAS tables
    ADASInfo = input.ADASInfo;
else
    ADASInfo = load('/home/verhaegh/DeuteriumADAS'); %load ADAS data set
end

%Uncertainty accounting in atomic rate coefficents
if isfield(input, 'CoeffUncertainty')
   PR = input.CoeffUncertainty.PR;
   PE = input.CoeffUncertainty.PE;
   PIR = input.CoeffUncertainty.PIR;
   PRR = input.CoeffUncertainty.PRR;
   PCX = input.CoeffUncertainty.PCX;
   PEP = input.CoeffUncertainty.PEP;
   PRP = input.CoeffUncertainty.PRP;
else
    PR = ones(numel(IR),1);
    PE = ones(numel(IR),1);
    PRR = ones(numel(IR),1);
    PIR = ones(numel(IR),1);
    PEP = ones(numel(IR),1);
    PRP = ones(numel(IR),1);
    PCX = ones(numel(IR),1);
end

%% Initialise matrix

TeEV = zeros(numel(IR),1)+NaN;
TeRV = zeros(numel(IR),1)+NaN;
RecV = zeros(numel(IR),1)+NaN;
IonV = zeros(numel(IR),1)+NaN;
PIonV = zeros(numel(IR),1)+NaN;
PRecV = zeros(numel(IR),1)+NaN;
CCV = zeros(numel(IR),1)+NaN;

%% Perform calculation
parfor i=1:numel(IE)
    
    if ~isnan(Den(i))&& ~isnan(none(i)) && ~isnan(IR(i)) && ~isnan(IE(i)) && ~isnan(DL(i)) %parameter check, if parameters not valid return NaN
        
    [ACD, SCD, PECrec, PECexc, PLT, PRB, CCD] = TECPEC(n1-2, Den(i)*ones(numel(Te),1), Te, ADASInfo, 1:7); %Interpolate ADAS parameters
    
    %modify rates depending on their uncertainty input
    PECrec = PR(i).*PECrec;
    PECexc = PECexc.*PE(i);
    SCD = SCD.*PRR(i);
    ACD = ACD.*PIR(i);
    PLT = PLT.*PEP(i);
    PRB = PRB.*PRP(i);
    CCD = CCD.*PCX(i);
    
    %perform intersections to get plasma parameters
    [~, TeR] = DSS_Aux_intersections(DL(i).*(Den(i)^2).*PECrec(:),Te, [1 1]*IR(i), [-Inf Inf]);
    [~, Rec] = DSS_Aux_intersections(DL(i).*(Den(i)^2).*PECrec(:),DL(i)*(Den(i)^2).*PECrec(:).*(SCD./PECrec(:)),[1 1]*IR(i), [-Inf Inf]);
    [~, PRec] = DSS_Aux_intersections(DL(i).*(Den(i)^2).*PECrec(:),DL(i).*(Den(i)^2).*PRB, [1 1]*IR(i), [-Inf Inf]); 
    [~, TeE] = DSS_Aux_intersections(DL(i).*((Den(i)^2).*none(i)).*PECexc(:), Te, [1 1]*IE(i), [-Inf Inf]);
    [~, Ion] = DSS_Aux_intersections(DL(i).*((Den(i)^2).*none(i)).*PECexc(:), DL(i).*((Den(i)^2).*none(i)).*PECexc(:).*(ACD./PECexc(:)), [1 1]*IE(i), [-Inf Inf]);    
    [~, PIon] = DSS_Aux_intersections(DL(i).*((Den(i)^2).*none(i)).*PECexc(:), DL(i).*((Den(i)^2).*none(i)).*PLT, [1 1]*IE(i), [-Inf Inf]);
    [~, CC] = DSS_Aux_intersections(DL(i).*((Den(i)^2).*none(i)).*PECexc(:), DL(i).*((Den(i)^2).*none(i)).*CCD, [1 1]*IE(i), [-Inf Inf]);

    
    if isempty(TeR)
        TeR = NaN;
        Rec = NaN;
        PRec = NaN;
    end
    
    if isempty(TeE)
        TeE = NaN;
        Ion = NaN;
        PIon = NaN;
    end
    
    if isempty(CC)
        CC = NaN;
    end
    
    if numel(TeE)>1
        TeE = NaN;
        Ion = NaN;
        PIon = NaN;
    end
    
    if numel(CC)>1
        CC = NaN;
    end
    
    if numel(TeR)>1
        TeR = NaN;
        Rec = NaN;
        PRec = NaN;
    end
    
    %% perform checks on output
    else
        TeE= NaN; %put NAN if values are invalid
        TeR= NaN;
        Rec = NaN;
        Ion = NaN;
        PRec = NaN;
        PIon = NaN;
        CC = NaN;
    end
    
    %% Put in matrix
    
    TeRV(i) = TeR;
    TeEV(i) = TeE;
    RecV(i) = Rec;
    IonV(i) = Ion;
    PRecV(i) = PRec;
    PIonV(i) = PIon;
    CCV(i) = CC;
end

%% Output calculation into output structure

output.Rec = reshape(RecV, size(IR));
output.Ion = reshape(IonV, size(IR));
output.TeR = reshape(TeRV, size(IR));
output.TeE = reshape(TeEV, size(IR));
output.PRec = reshape(PRecV, size(IR));
output.PIon = reshape(PIonV, size(IR));
output.CC = reshape(CCV, size(IR));
