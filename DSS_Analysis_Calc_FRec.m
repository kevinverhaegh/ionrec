function output = DSS_Analysis_Calc_FRec(input)
%Function to determine FRec 

%% Analyse input, get necessary info from input structure

n1 = input.n1; %First balmer line index
n2 = input.n2; %second Balmer line index
if n1>n2
    error('n2 should correspond to the higher Balmer line -> please reverse')
end
n1Int = input.n1Int; %First Balmer line intensity
n2Int = input.n2Int; %second Balmer line intensity
none = input.none; %neutral fraction
Den = input.Den; %electron density

%% Get optional info from input structure
if isfield(input, 'Te') %initial Te vector
    Te = input.Te;
else
    Te = 0.5:0.01:100; 
end

if isfield(input, 'ADASInfo') %ADAS tables
    ADASInfo = input.ADASInfo;
else
    ADASInfo = load('/home/verhaegh/DeuteriumADAS'); %load ADAS data set
end

if isfield(input, 'TeRestrict') %Te restrict implies whether to cut-off the used Te vector to elliminate 
    TeRestrict = input.TeRestrict;
else
    TeRestrict = 1;
end

if isfield(input, 'ExtremaSol') %maximum difference for intersect values to be lumped
    ExtremaSol = input.ExtremaSol;
else
    ExtremaSol = 0;
end

if isfield(input, 'Dupl') %maximum number of curve intersections
    Dupl = input.Dupl;
else
    Dupl = 1;
end 

if isfield(input, 'CoeffUncertainty') %Uncertainties in the atomic data
   PR1 = input.CoeffUncertainty.PR1;
   PE1 = input.CoeffUncertainty.PE1;
   PR2 = input.CoeffUncertainty.PR2;
   PE2 = input.CoeffUncertainty.PE2;
else
    PR1 = ones(numel(n1Int),1);
    PE1 = ones(numel(n1Int),1);
    PE2 = ones(numel(n1Int),1);
    PR2 = ones(numel(n1Int),1);
end

%% Initialise cell for storing result

TeC = cell(numel(n1Int),1);
FRec1C = cell(numel(n1Int),1);
FRec2C = cell(numel(n1Int),1);
ExtremaLR = zeros(numel(n1Int),2);
ExtremaFR1 = zeros(numel(n1Int),2);
ExtremaFR2 = zeros(numel(n1Int),2);

TeO = Te;

%% Perform calculation
parfor i=1:numel(n1Int)
    if ~isnan(Den(i))&& ~isnan(none(i)) && ~isnan(n1Int(i)) && ~isnan(n2Int(i)) %parameter check, if parameters not valid return NaN
        
    Te=TeO;
    [~, ~, PECrec, PECexc] = TECPEC([n1 n2]-2, Den(i)*ones(size(Te)), Te, ADASInfo,3:4); %Load ADAS coefficients
    PECrec(1,:) = PR1(i).*PECrec(1,:);
    PECrec(2,:) = PR2(i).*PECrec(2,:);
    PECexc(1,:) = PE1(i).*PECexc(1,:);
    PECexc(2,:) = PE2(i).*PECexc(2,:);
    
    LineRat = squeeze(PECrec(2,:,:) + none(i).*PECexc(2,:,:))./squeeze(PECrec(1,:,:) + none(i).*PECexc(1,:,:)); %model line ratio
    
    if TeRestrict %restrict Te based on the minima and maxima of the Balmer line ratio - if requested
        [~,TeE] = min(LineRat);
        [~,TeB] = max(LineRat);
        Te=Te(TeB:TeE);
        LineRat=LineRat(TeB:TeE);
        PECrec = PECrec(:,TeB:TeE);
        PECexc = PECexc(:,TeB:TeE);
        FRec1Min = min(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:)));
        FRec1Max = max(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:)));
        FRec2Min = min(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:)));
        FRec2Max = max(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:)));        
    end

    [~, TeS] = DSS_Aux_intersections(LineRat, Te, [1 1].*(n2Int(i)./n1Int(i)), [-inf inf]); %obtain Te at which the modelled line ratio matches the measured line ratio
    [~, FRec1S] = DSS_Aux_intersections(LineRat, PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:)), [1 1].*(n2Int(i)./n1Int(i)), [-inf inf]);
    [~, FRec2S] = DSS_Aux_intersections(LineRat, PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:)), [1 1].*(n2Int(i)./n1Int(i)), [-inf inf]);

    %deprecated
    if ExtremaSol
        if numel(TeS)>1
            %take extram solutions
            if FRec1S<0.5
                FRec1S = mean(FRec1S);
                FRec2S = mean(FRec2S);
                TeS = mean(TeS);
            else
                FRec1S = mean(FRec1S);
                FRec2S = mean(FRec2S);
                TeS = mean(TeS);
            end
        end
    end

%% Perform checks on output

    if isempty(FRec1S) %if no solution
        if (n2Int(i)/n1Int(i))>max(LineRat) %checked if the measured line ratio is higher than the maximum line ratio
            FRec1S = max(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:)))  + rand()*(1-max(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:)))); %if so, use random value between maximum determinable FRec and 1
            FRec2S = max(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:)))  + rand()*(1-max(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:))));
            TeS = min(Te);
        elseif (n2Int(i)/n1Int(i))<min(LineRat) %checked if the measured line ratio is lower than the minimum line ratio
            FRec1S = rand()*min(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:))); %if so, use random value between maximum determinable FRec and 1
            FRec2S = rand()*min(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:)));
            TeS = max(Te);
         end
     else
%          EF1 = ((n2Int(i)/n1Int(i)) - min(LineRat))/((FRec1Max - FRec1Min).*((n2Int(i)/n1Int(i)) - min(LineRat)) + FRec1Min.*(max(LineRat) - min(LineRat)));
%          FRec1S = FRec1S.*EF1;
%          EF2 = ((n2Int(i)/n1Int(i)) - min(LineRat))/((FRec2Max - FRec2Min).*((n2Int(i)/n1Int(i)) - min(LineRat)) + FRec2Min.*(max(LineRat) - min(LineRat)));
%          FRec2S = FRec2S.*EF2;
    end
    
    ExtremaLR(i,:) = [min(LineRat) max(LineRat)]; 
    ExtremaFR1(i,:) = [min(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:))) max(PECrec(1,:)./(PECrec(1,:) + none(i).*PECexc(1,:)))]; 
    ExtremaFR2(i,:) = [min(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:))) max(PECrec(2,:)./(PECrec(2,:) + none(i).*PECexc(2,:)))]; 
    

    else
        TeS = NaN; %put NAN if values are invalid
        FRec1S = NaN;
        FRec2S = NaN;
    end
    %% Put in matrix

    TeC{i} = TeS;
    FRec1C{i} = FRec1S;
    FRec2C{i} = FRec2S;
end

%% transforming cell into matrix


TeV = zeros(numel(n1Int),Dupl)+NaN;
FRec1V = zeros(numel(n1Int),Dupl)+NaN;
FRec2V = zeros(numel(n1Int),Dupl)+NaN;

for i=1:numel(n1Int)
    if ~isempty(TeC{i})
        if numel(TeC{i})<=Dupl
            TeV(i,1:numel(TeC{i})) = TeC{i};
            FRec1V(i,1:numel(TeC{i})) = FRec1C{i};
            FRec2V(i,1:numel(TeC{i})) = FRec2C{i};
        end
    end
end

%% Output calculation into output structure

output.FRec1 =reshape(FRec1V,size(n1Int));
output.FRec2 = reshape(FRec2V,size(n1Int));
output.Te = reshape(TeV,size(n1Int));
output.Extrema.ExtremaLR = ExtremaLR;
output.Extrema.ExtremaFR1 = ExtremaFR1;
output.Extrema.ExtremaFR2 = ExtremaFR2;