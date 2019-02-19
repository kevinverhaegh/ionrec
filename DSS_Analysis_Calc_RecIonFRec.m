function output = DSS_Analysis_Calc_RecIonFRec(input)
%Kevin Verhaegh, University of York / SPC-EPFL / CCFE - UK
%kevin.verhaegh@ukaea.uk

%Code to estimate ionisation/recombination rates from Balmer line
%emissivities in a Monte-Carlo probablistic fashion. See pre-print X.

%This code requires a specific ADAS structure, which for legal reasons
%cannot be provided with the code. Please either ask the author or note the
%comments below to obtain a matlab structure for this ADAS format

%This code requires the Matlab image processing, statistics and parallel
%computation toolboxes.

%% Disable auto parpool creation
 ps = parallel.Settings;
  ps.Pool.AutoCreate = false;
 delete(gcp('nocreate')); %close off parpool

%% Analyse input, get necessary info from input structure
%username=char(java.lang.System.getProperty('user.name'));

floc = input.floc; %file save location
disp('Analysing input....')

%floc = ['/NoTivoli/',username,'/DSS_Analysis/', num2str(input.shot), '_', num2str(input.n1), num2str(input.n2), '.mat' ];
save(floc, 'input', '-v7.3');

n1 = input.n1; %First balmer line index
n2 = input.n2; %second Balmer line index
if n1>n2
    error('n2 should correspond to the higher Balmer line -> please reverse')
end
n1Int = input.n1Int; %First Balmer line intensity
n1IntErr = input.n1IntErr; %1-sigma error First Balmer line intensity (random)
n2Int = input.n2Int; %second Balmer line intensity
RelErr = input.RelErr;
none = input.none; %neutral fraction
noneL = input.noneL; %lower limit neutral fraction
noneH = input.noneH; %upper limit neutral fraction
Den = input.Den; %electron density
DenErr = input.DenErr; %random error electron density
DL = input.DL; %Divertor leg length
DLL = input.DLL; %lower estimate divertor leg length
DLH = input.DLH; %upper estimate divertor leg length
R = input.R; %radius 
Z = input.Z; %z-position
Time = input.Time; %time vector
shot = input.shot; %shot number

%% Get optional info from input structure

%number of random parameters
if isfield(input, 'Iter')
    Iter = input.Iter;
else
    Iter = 5000;
end

%ADAS
if isfield(input, 'ADASInfo')
    ADASInfo = input.ADASInfo;
else
    error('An overview of ADAS data needs to be provided with the analysis - in a particular format (see comment in code or ask author)');
    %ADASInfo = load('/home/verhaegh/DeuteriumADAS');
    %ADAS structure:
    %ADASInfo.AdasNe %ADAS density grid
    %ADASInfo.AdasTe %ADAS temperature grid
    %ADASInfo.SCDHydrogen %Effective ionisation rates [density x temperature]
    %ADASInfo.ACDHYdrogen %Effective recombination rates [density x
    %temperature]
    %ADASInfo.BalmerRecombination %Recombination PEC [10 x temperature x density]
    %1-10 depicts different Balmer lines (n=3,4,5,6,...)
    %ADASInfo.BalmerExcitation %Excitation PEC [10 x temperature x density]
    %1-10 depicts different Balmer lines (n=3,4,5,6,...)
    %ADASInfo.PLTHydrogen %Hydrogen radiative power losses - ADAS PLT
    %coefficient [density x temperature]
    %ADASInfo.PRBHydrogen %Hydrogen radiative power losses due to
    %recombination - ADAS PRB [density x temperature]
    %Data should be provided from Open-ADAS 2012 dataset and provided in
    %the above structure.
end

%Minimum density
if isfield(input, 'DenMin')
    DenMin = input.DenMin;
else
    DenMin = 5e18;
end

%Correction algorithm for Frec
if isfield(input, 'FRecForceTrend')
    DoFRecForceTrend = input.FRecForceTrend;
else
    DoFRecForceTrend = 1;
end

%Te vector used for FRec determination
if isfield(input, 'TeRestrict')
    TeRestrict = input.TeRestrict;
else
    TeRestrict = 1;
end

%Deprecated
if isfield(input, 'ExtremaSol')
    ExtremaSol = input.ExtremaSol;
else
    ExtremaSol = 0;
end

%Quantile probability calculation (Deprecated)
if isfield(input, 'Quant')
    Quant = input.Quant;
else
    Quant = linspace(0, 1, 100);
end

%Parallel execution yes/no
if isfield(input, 'DoParFor')
    if input.DoParFor == 1
        parpool()
    end
end

%Correction algorithm for increase in TeE
if ~isfield(input, 'TeE_filter')
    input.TeE_filter = 1;
end

%Time at which the correction algorithm starts
if ~isfield(input, 'TeE_filter')
    if ~isnfield(input, 'TeE_filterTime')
        input.TeE_filterTime = 1;
    end
end

%Inclusion of uncertainties in the ADAS coefficients
if isfield(input, 'CoeffUncertainty')
    DoCoeffUncertainty = 1;
    CoeffUncertainty = input.CoeffUncertainty;
else
    DoCoeffUncertainty = 0;
    CoeffUncertainty = 0;
end

%Log-uniform distribution (as opposed to the default uniform) for no/ne
if isfield(input,'none_loguniform')
    none_loguniform = input.none_loguniform;
else
    none_loguniform = 0;
end

%% Initialise random values

disp('Intialising matrices and random values....')

%make empty cells

n1IntMC = zeros([size(n1Int),Iter])+NaN;
n2IntMC = zeros([size(n1Int),Iter])+NaN;
noneMC = zeros([size(n1Int),Iter])+NaN;
DenMC = zeros([size(n1Int), Iter]) + NaN;
DLMC = zeros([size(n1Int),Iter])+NaN;

%initialise random vectors

noneRand = rand(Iter,1); %random vector none
randn_n1 = randn(Iter,1);    
randn_LRMC = randn(Iter, 1);
randn_DenMC = randn(Iter, 1);
randn_DL = randn(Iter, 1);

%random vector for rate uncertainty
PR1 = abs(1 + CoeffUncertainty.*rand(Iter,1));
PE1 = abs(1 + CoeffUncertainty.*rand(Iter,1));
PR2 = abs(1 + CoeffUncertainty.*rand(Iter,1));
PE2 =abs(1 + CoeffUncertainty.*rand(Iter,1));
PIR = abs(1 + CoeffUncertainty.*rand(Iter,1));
PRR = abs(1 + CoeffUncertainty.*rand(Iter,1));
PCX = abs(1 + CoeffUncertainty.*rand(Iter,1));
PEP = abs(1 + CoeffUncertainty.*rand(Iter,1));
PRP = abs(1 + CoeffUncertainty.*rand(Iter,1));

%adjust randn_DenMC
%% utilize random vectors to assign the MC values for the density using rejection samplign to account for minimum density
for j=1:numel(n1Int(1,:)) 
    for i=1:numel(n1Int(:,1))
        if ~isnan(Den(i,j))&& ~isnan(none(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(DL(i,j)) %parameter check, if parameters not valid return NaN

            DenMC(i,j,:) = Den(i,j) + DenErr(i,j).*randn_DenMC(:);
            NOK = 1;
            while NOK==1
               
                Slice = squeeze(DenMC(i,j,:));
                V = find(Slice<DenMin);
                if ~isempty(V)
                Vr = randn(numel(V),1);
                Vn = Den(i,j) + DenErr(i,j).*Vr; 
                randn_DenMC(V) = Vr;
                DenMC(i,j,(DenMC(i,j,:)<DenMin)) = Vn;
                else
                    NOK=0;
                end
            end

        end
    end
end
DenMC = zeros([size(n1Int),Iter])+NaN;
%% utilize random vectors to assign all the other MC values
for j=1:numel(n1Int(1,:)) 
    for i=1:numel(n1Int(:,1))
        if ~isnan(Den(i,j))&& ~isnan(none(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(DL(i,j)) %parameter check, if parameters not valid return NaN

            n1IntMC(i,j,:) = n1Int(i,j)*(1 + randn_n1(:).*(n1IntErr(i,j)/n1Int(i,j)));
            LRMC = squeeze((n2Int(i,j)./n1Int(i,j))).*(1+RelErr.*randn_LRMC(:));
            n2IntMC(i,j,:) = squeeze(n1IntMC(i,j,:)).*LRMC;
            %n2IntMC(i,j,:)  = n2Int(i,j).*(1 + n2Int(i,j)./n1Int(i,j))randn(Iter,1))n2Int(i,j)*(1 + randn(Iter,1).*(n2IntErr(i,j)/n2Int(i,j)));
            if ~none_loguniform
            noneMC(i,j,:)  = (noneH(i,j) - noneL(i,j))*noneRand + noneL(i,j);
            else
            noneMC(i,j,:) = DSS_Aux_LogUniform_Sample(noneRand, noneL(i,j), noneH(i,j));     
            end
            DenMC(i,j,:) = Den(i,j) + DenErr(i,j).*randn_DenMC(:);
            
            randV = randn_DL(:);
            randV(randV<0) = (DL(i,j) - DLL(i,j)).*randV(randV<0);
            randV(randV>0) = (DLH(i,j) - DL(i,j)).*randV(randV>0);
            DLMC(i,j,:) = DL(i,j) + randV;
            

        end
    end
end



%% Make empty matrices for storing calculation results

%main results
FRec1 = 0.*n1Int + NaN;
FRec2 = 0.*n1Int + NaN;
Te = 0.*n1Int + NaN;
Rec1 = 0.*n1Int + NaN;
Ion1 = 0.*n1Int + NaN;
TeR1 = 0.*n1Int + NaN;
TeE1 = 0.*n1Int + NaN;
PIon1 = 0.*n1Int + NaN;
PRec1 = 0.*n1Int + NaN;
CX1 = 0.*n1Int + NaN;
Rec2 = 0.*n1Int + NaN;
Ion2 = 0.*n1Int + NaN;
TeR2 = 0.*n1Int + NaN;
TeE2 = 0.*n1Int + NaN;
PIon2 = 0.*n1Int + NaN;
PRec2 = 0.*n1Int + NaN;
CX2 = 0.*n1Int + NaN;

%monte carlo solutions
FRec1MC = 0.*n1IntMC + NaN;
FRec2MC = 0.*n1IntMC + NaN;
TeMC = 0.*n1IntMC + NaN;
RecMC1 = 0.*n1IntMC + NaN;
IonMC1 = 0.*n1IntMC + NaN;
TeRMC1 = 0.*n1IntMC + NaN;
TeEMC1 = 0.*n1IntMC + NaN;
PIonMC1 = 0.*n1IntMC + NaN;
PRecMC1 = 0.*n1IntMC + NaN;
CXMC1 = 0.*n1IntMC + NaN;
RecMC2 = 0.*n1IntMC + NaN;
IonMC2 = 0.*n1IntMC + NaN;
TeRMC2 = 0.*n1IntMC + NaN;
TeEMC2 = 0.*n1IntMC + NaN;
PIonMC2 = 0.*n1IntMC + NaN;
PRecMC2 = 0.*n1IntMC + NaN;
CXMC2 = 0.*n1IntMC + NaN;

%% initialise settings for FRec calculation

%make general input settings
inputSet.n1 = n1;
inputSet.n2 = n2;
inputSet.ADASInfo = ADASInfo;
inputSet.Dupl = 1;
inputSet.TeRestrict = 1;
inputSet.ExtremaSol = 0;

%store general input settings in inputVal (main result) and inputMC (monte
%carlo result)
inputVal = inputSet;
inputMC = inputSet;

dimen=size(n1Int);

%% Start FRec calculation loop
disp('Starting FRec calculation....')
looperT = 0; %get index of values to fill in for estimating the duration of this step
for i=1:dimen(1)
    for j=1:dimen(2)
        if ~isnan(Den(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(none(i,j)) 
            looperT=looperT+1;
        end
    end
end
tic
looper=0;
for i=1:dimen(1)
   for j=1:dimen(2)
    
        %perform parameter check
        if ~isnan(Den(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(none(i,j)) 
            looper=looper+1;
            %set-up inputVal (for determining FRec)
            inputVal.n1Int = n1Int(i,j);
            inputVal.n2Int = n2Int(i,j);
            inputVal.none = none(i,j);
            inputVal.Den = Den(i,j);
            inputVal.TeRestrict = TeRestrict;
            inputVal.ExtremaSol = ExtremaSol;
            inputVal.DoParFor = input.DoParFor;
            
            outputVal = DSS_Analysis_Calc_FRec(inputVal);

            %set-up inputMC
            inputMC.n1Int = n1IntMC(i,j,:);
            inputMC.n2Int = n2IntMC(i,j,:);
            inputMC.none = noneMC(i,j,:);
            inputMC.Den = DenMC(i,j,:);
            inputMC.TeRestrict = TeRestrict;
            inputMC.ExtremaSol = ExtremaSol;
            inputMC.DoParFor = input.DoParFor;
            
            if DoCoeffUncertainty
                inputMC.CoeffUncertainty.PR1 = PR1;
                inputMC.CoeffUncertainty.PE1 = PE1;
                inputMC.CoeffUncertainty.PR2 = PR2;
                inputMC.CoeffUncertainty.PE2 = PE2;
            end

            outputMC = DSS_Analysis_Calc_FRec(inputMC);

            %obtain data from outputVal and put it in a matrix
            
            FRec1(i,j) = outputVal.FRec1;
            FRec2(i,j) = outputVal.FRec2;
            Te(i,j) = outputVal.Te;
            
            %obtain data from outputMC and put it in a cell
            FRec1MC(i,j,:) = squeeze(outputMC.FRec1);
            FRec2MC(i,j,:) = squeeze(outputMC.FRec2);
            TeMC(i,j,:) = squeeze(outputMC.Te);
            
            
            tp = toc;
            disp(['Estimated time left: ', num2str((tp/looper)*looperT - tp)])
            
        end
        
    end
end
tp = toc;
disp(['FRec calculation finished - ', num2str(tp), ' seconds elapsed'])

%% Perform post-processing on FRec calculation to clean-up the FRec vaectors

%first force FRec trend - if requested (correction algorithm to ensure
%correct FRec values at high FRec.
if DoFRecForceTrend
    disp('Performing post-processing on FRec - enforcing increasing FRec trend....')
    %make backups of calculated data
    FRec1MC_o = FRec1MC;
    FRec2MC_o = FRec2MC;
    TeMC_o = TeMC;
    FRec1_o = FRec1;
    FRec2_o = FRec2;
    Te_o = Te;
    for i=1:numel(n1Int(:,1))
        for j=1:Iter
            FRec1MC(i,:,j) = FRecForceTrend(FRec1MC(i,:,j)); %apply post-process FRec
            FRec2MC(i,:,j) = FRecForceTrend(FRec2MC(i,:,j));
            TeMC(i,:,j) = FRecForceTrend(TeMC(i,:,j));
        end
        %deprecated
        FRec1(i,:) = FRecForceTrend(FRec1(i,:));
        FRec2(i,:) = FRecForceTrend(FRec2(i,:));
        Te(i,:) = FRecForceTrend(Te(i,:));
    end
    disp('FRec post-processing finished')
    
    %deprecated
 FRec1 = DSS_Aux_smooth2a(FRec1, 2,2); %smooth post-processed FRec
 FRec1 = imadjust(FRec1, [min(FRec1(:)) max(FRec1(:))], [min(FRec1_o(:)) max(FRec1_o(:))]); %adjust post-smoothed FRec to realign with previous minima/maxima
% 
 FRec2 = DSS_Aux_smooth2a(FRec2, 2,2);
 FRec2 = imadjust(FRec2, [min(FRec2(:)) max(FRec2(:))], [min(FRec2_o(:)) max(FRec2_o(:))]); %adjust post-smoothed FRec to realign with previous minima/maxima
% 
 for i=1:Iter %now apply smoothing of post-processed FRec on the Monte Carlo execution grid
     MiP = min(min(min(FRec1MC(:,:,i))));
     MaP = max(max(max(FRec1MC(:,:,i))));
     
     FRec1MC(:,:,i) = FRec1MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
     
     FRec1MC(:,:,i) =  DSS_Aux_smooth2a(squeeze(FRec1MC(:,:,i)), 2,2);   %Apply smoothing
     FRec1MC(:,:,i) = imadjust(squeeze(FRec1MC(:,:,i)), [min(min(min(FRec1MC(:,:,i)))) max(max(max(FRec1MC(:,:,i))))], [MiP MaP]); %ajust smoothed FRec to make it align with previous minima/maxima
     
     FRec1MC(:,:,i) = FRec1MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
     
     MiP = min(min(min(FRec2MC(:,:,i))));
     MaP = max(max(max(FRec2MC(:,:,i))));
     
     FRec2MC(:,:,i) = FRec2MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
     
     FRec2MC(:,:,i) =  DSS_Aux_smooth2a(squeeze(FRec2MC(:,:,i)), 2,2);  %Apply smoothing 
     FRec2MC(:,:,i) = imadjust(squeeze(FRec2MC(:,:,i)), [min(min(min(FRec2MC(:,:,i)))) max(max(max(FRec2MC(:,:,i))))], [MiP MaP]); %ajust smoothed FRec to make it align with previous minima/maxima
    
     FRec2MC(:,:,i) = FRec2MC(:,:,i).*(n1Int./n1Int).*(n2Int./n2Int); %Get FFRec filtering NaN for both intensity matrices
   
 end
 
    
end
% 

%% Perform recombination/ionisation rate calculation

%make general input settings
clear inputSet
inputSet.ADASInfo = ADASInfo;

%store general input settings in inputVal (main result) and inputMC (monte
%carlo result)
inputVal = inputSet;
inputMC = inputSet;
disp('Starting recombination / ionisation rate calculations for both Balmer lines...')
tic
looper=0;
for i=1:numel(n1Int(:,1))
    for j=1:numel(n1Int(1,:))
          
        if ~isnan(Den(i,j)) && ~isnan(n1Int(i,j)) && ~isnan(n2Int(i,j)) && ~isnan(none(i,j)) 
            looper=looper+1;
        %set-up inputVal for first Balmer line index
        inputVal.n1 = n1;
        inputVal.IR = n1Int(i,j).*FRec1(i,j);
        inputVal.IE = n1Int(i,j).*(1-FRec1(i,j));
        inputVal.none = none(i,j);
        inputVal.Den = Den(i,j);
        inputVal.DL = DL(i,j);
        inputVal.ADASInfo = ADASInfo;
        
        outputValn1 = DSS_Analysis_Calc_RecIon(inputVal);
        
        %set-up inputMC for first Balmer line index
        inputMC.n1 = n1;
        inputMC.IR = n1IntMC(i,j,:).*FRec1MC(i,j,:);
        inputMC.IE = n1IntMC(i,j,:).*(1-FRec1MC(i,j,:));
        inputMC.none = noneMC(i,j,:);
        inputMC.Den = DenMC(i,j,:);
        inputMC.DL = DLMC(i,j,:);
        inputMC.ADASInfo = ADASInfo;
        
        if DoCoeffUncertainty %Uncertainty parameters in atomic rate coefficients
        inputMC.CoeffUncertainty.PR = PR1;
        inputMC.CoeffUncertainty.PE = PE1;
        inputMC.CoeffUncertainty.PIR = PIR;
        inputMC.CoeffUncertainty.PRR = PRR;
        inputMC.CoeffUncertainty.PCX = PCX;
        inputMC.CoeffUncertainty.PEP = PEP;
        inputMC.CoeffUncertainty.PRP = PRP;
        end
        
        outputMCn1 = DSS_Analysis_Calc_RecIon(inputMC); %Perform ionisation/recombination calculation
           
        %obtain data from outputVal and put it in a matrix
        Rec1(i,j) = outputValn1.Rec;
        Ion1(i,j) = outputValn1.Ion;
        TeE1(i,j) = outputValn1.TeE;
        TeR1(i,j) = outputValn1.TeR;
        PIon1(i,j) = outputValn1.PIon;
        PRec1(i,j) = outputValn1.PRec;
        CX1(i,j) = outputValn1.CC;
    
        %obtain data from outputMC and put it in a cell
        RecMC1(i,j,:) = outputMCn1.Rec;
        IonMC1(i,j,:) = outputMCn1.Ion;
        TeEMC1(i,j,:) = outputMCn1.TeE;
        TeRMC1(i,j,:) = outputMCn1.TeR;
        PIonMC1(i,j,:) = outputMCn1.PIon;
        PRecMC1(i,j,:) = outputMCn1.PRec;
        CXMC1(i,j,:) = outputMCn1.CC;
      
        %set-up inputVal for second Balmer line index
        inputVal.n1 = n2;         
        inputVal.IR = n2Int(i,j).*FRec2(i,j);
        inputVal.IE = n2Int(i,j).*(1-FRec2(i,j));
        inputVal.none = none(i,j);
        inputVal.Den = Den(i,j);
        inputVal.DL = DL(i,j);
        inputVal.ADASInfo = ADASInfo;
        
        outputValn2 = DSS_Analysis_Calc_RecIon(inputVal);
        
        %set-up inputMC for second Balmer line index
        inputMC.n1 = n2;     
        inputMC.IR = n2IntMC(i,j,:).*FRec2MC(i,j,:);
        inputMC.IE = n2IntMC(i,j,:).*(1-FRec2MC(i,j,:));
        inputMC.none = noneMC(i,j,:);
        inputMC.Den = DenMC(i,j,:);
        inputMC.DL = DLMC(i,j,:);
        inputMC.ADASInfo = ADASInfo;
        
        if DoCoeffUncertainty %Uncertainties atomic rate coefficients
        inputMC.CoeffUncertainty.PR = PR2;
        inputMC.CoeffUncertainty.PE = PE2;
        inputMC.CoeffUncertainty.PIR = PIR;
        inputMC.CoeffUncertainty.PRR = PRR;
        inputMC.CoeffUncertainty.PCX = PCX;
        inputMC.CoeffUncertainty.PEP = PEP;
        inputMC.CoeffUncertainty.PRP = PRP;
        end
        
        outputMCn2 = DSS_Analysis_Calc_RecIon(inputMC); %Get ionisation/recombination rates
           
        %obtain data from outputVal and put it in a matrix
        Rec2(i,j) = outputValn2.Rec;
        Ion2(i,j) = outputValn2.Ion;
        TeE2(i,j) = outputValn2.TeE;
        TeR2(i,j) = outputValn2.TeR;
        PIon2(i,j) = outputValn2.PIon;
        PRec2(i,j) = outputValn2.PRec;
        CX2(i,j) = outputValn2.CC;
    
        %obtain data from outputMC and put it in a cell
        RecMC2(i,j,:) = outputMCn2.Rec;
        IonMC2(i,j,:) = outputMCn2.Ion;
        TeEMC2(i,j,:) = outputMCn2.TeE;
        TeRMC2(i,j,:) = outputMCn2.TeR;  
        PIonMC2(i,j,:) = outputMCn2.PIon;
        PRecMC2(i,j,:) = outputMCn2.PRec;
        CXMC2(i,j,:) = outputMCn2.CC;
        
            tp = toc;
            disp(['Estimated time left: ', num2str((tp/looper)*looperT - tp)])       
        
        end
        
    end
end

tp = toc;
disp(['Rec/Ion calculation finished - ', num2str(tp), ' seconds elapsed'])
  
%% Wrap up results in output structure 

%input matrices and standard parameters
output.input = input;
output.input.Iter = Iter;
output.input.ADASInfo = ADASInfo;
output.input.DenMin = DenMin;
output.input.FRecForceTrend = DoFRecForceTrend;
output.input.Quant = Quant;

%MC input
output.inputMC.n1Int = n1IntMC;
output.inputMC.n2Int = n2IntMC;
output.inputMC.none = noneMC;
output.inputMC.DL = DLMC;
output.inputMC.Den = DenMC;

%main results - FRec
output.Result.FRec1 = FRec1;
output.Result.FRec2 = FRec2;
output.Result.Te = Te;

%main results - Rec
output.Result.Rec1 = Rec1;
output.Result.PRec1 = PRec1;
output.Result.TeR1 = TeR1;

output.Result.Rec2 = Rec2;
output.Result.PRec2 = PRec2;
output.Result.TeR2 = TeR2;

%main results - Ion
output.Result.Ion1 = Ion1;
output.Result.PIon1 = PIon1;
output.Result.TeE1 = TeE1;
output.Result.CX1 = CX1;

output.Result.Ion2 = Ion2;
output.Result.PIon2 = PIon2;
output.Result.TeE2 = TeE2;
output.Result.CX2 = CX2;

%main results - FRec - MC
output.ResultMC.FRec1MC = FRec1MC;
output.ResultMC.FRec2MC = FRec2MC;
output.ResultMC.TeMC = TeMC;

%main results - Ion - MC
output.ResultMC.Ion1MC = IonMC1;
output.ResultMC.PIon1MC = PIonMC1;
output.ResultMC.TeE1MC = TeEMC1;
output.ResultMC.Ion2MC = IonMC2;
output.ResultMC.PIon2MC = PIonMC2;
output.ResultMC.TeE2MC = TeEMC2;
output.ResultMC.CX1MC = CXMC1;
output.ResultMC.CX2MC = CXMC2;

%main results - Rec - MC
output.ResultMC.Rec1MC = RecMC1;
output.ResultMC.PRec1MC = PRecMC1;
output.ResultMC.TeR1MC = TeRMC1;
output.ResultMC.Rec2MC = RecMC2;
output.ResultMC.PRec2MC = PRecMC2;
output.ResultMC.TeR2MC = TeRMC2;

%try to get bolo / lp analysis  
% try
%     if shot < 80000
% output.BOLO = GetBoloPW(output.input.shot, 0.6);
% output.LP = GetFluxSP2(output.input.shot);
%     end
% catch
% end

%save parameters locally
%clearvars -EXCEPT output username
%floc = ['/NoTivoli/', username, '/DSS_Analysis/', num2str(output.input.shot), '_', num2str(output.input.n1), num2str(output.input.n2), '.mat' ];
save(floc, 'output', '-v7.3');
%jheapcl %java heap cleanup
output = DSS_Analysis_Unpack(floc, output.input.TeE_filter, 1);
%java heap cleanup
try 
    jheapcl
catch
end

function FRecC = FRecForceTrend(FRec)
%When calculating FRec, during a density ramp experiment, FRec should be increasing, resulting in a hyperbolic tangent curve as function of time.
%However, due to the issues with converging to a unique FRec number at low/high FRec, this sometimes results in a 'S-like' behaviour.

%This function tends to rectify this behaviour

Window=0.1;

FRecC = FRec;

for i=1:numel(FRec(:,1))
    %for each wavelength
    Validity = ~isnan(FRec(i,:));
    if sum(Validity)>1
        %FRecTime = smooth(FRec(i,:),Window,'loess').*Validity';
        FRecTime = FRec(i,:);
        FRecTime(~Validity)=NaN;

        %get min/max of FrecTime
        [~,imin]=min(FRecTime);
        [~,imax] = max(FRecTime);

        indxmin = find(Validity);
        indxmin=indxmin(indxmin<imin);
        indxmax = find(Validity);
        indxmax=indxmax(indxmax>imax);

        FRecC(i,indxmin) = FRec(i,imin);
        FRecC(i,indxmax) = FRec(i,imax);
    
    end
    
end
