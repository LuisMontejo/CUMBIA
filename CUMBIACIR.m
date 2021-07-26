
%==========================================================================================

%                                    CUMBIACIR v0.2                 

%           SECTION AND MEMBER RESPONSE OF RC MEMBERS OF CIRCULAR SECTION

%                     LUIS A. MONTEJO (luis.montejo@upr.edu)

%              uptades available at https://github.com/LuisMontejo/CUMBIA

% cite as: Montejo, L. A., & Kowalsky, M. J. (2007). CUMBIA—Set of codes for the analysis of 
% reinforced concrete members. CFL technical rep. no. IS-07, 1.

%                            last updated: 07/26/2021

%==========================================================================================

clc;clear;close all;format long

% input data:

name = 'CUMBIACIREX' ;     %identifies actual work, the output file will be name.xls

interaction = 'y';      % if you want to also perform an axial load - moment interaction
                         % analysis type 'y', otherwise type 'n

% section properties:

D       = 1000;                          % section diameter (mm)
clb     = 50;                             % cover to longitudinal bars (mm)

% member properties

L             = 3000;                    % member clear length (mm)
bending       = 'single';                % single or double
ductilitymode = 'biaxial';              % biaxial or uniaxial

% reinforcement details:

nbl     = 22;             % number of longitudinal bars
Dbl     = 25;             % long. bar diameter (mm)   
Dh      = 9;              % diameter of transverse reinf. (mm)
type    = 'spirals';        % 'spirals' or 'hoops'*
s       = 120;               % spacing of transverse steel (mm)*

% aplieed loads:

P      =  400;                  % axial load kN (-) tension (+)compression

% material models (input the 'name' of the file with the stress-strain relationship
% to use the default models: Mander model for confined or unconfined  concrete type 'mc' or 'mu'.
% For lightweight confined concrete type 'mclw' 
% King model for the steel 'ks', Raynor model for steel 'ra':

confined   = 'mc';
unconfined = 'mu';
rebar      = 'ra';

% material properties 

fpc     = 35;                % concrete compressive strength (MPa)
Ec      = 0;                 % concrete modulus of elasticity (MPa) or
                             % input 0 for automatic calculation using
                             % 5000(fpc)^0.5
eco     = 0.002;            % unconfined strain (usually 0.002 for normal weight or 0.004 for lightweight)*
esm     = 0.11;              % max transv. steel strain (usually ~0.10-0.15)*
espall  = 0.0064;            % max uncon. conc. strain (usually 0.0064)

fy      = 460;               % long steel yielding stress (MPa)
fyh     = 400;               % transverse steel yielding stress (MPa)
Es      = 200000;            % steel modulus of elasticity
fsu     = 620;               % long steel max stress (MPa)*
esh     = 0.008;             % long steel strain for strain hardening (usually 0.008)*
esu     = 0.12;              % long. steel maximum strain (usually ~0.10-0.15)*

Ey     =  350;                % slope of the yield plateau (MPa)
C1     =  3.5;                % defines strain hardening curve in the Raynor model [2-6]

% *this information is used only if the default material models are selected

% strain limits for yield surface (interaction diagram);

csid = 0.004;  % concrete
ssid = 0.015;  % steel

% Deformation Limit States:

ecser = 0.004;      esser = 0.015;   % concrete (ecser) and steel (esser) serviceability strain
ecdam = 0.018;     esdam = 0.060;   % concrete (ecser) and steel (esser) damage control strain
                                   % (to use the 2/3 of the ultimate concrete strain just tipe 'twth'

% temperature information (in case of freezing conditions)
temp  = 40;           % temperature of the specimen in celsius
kLsp  = 0.022;        % constant to calculate Lsp = kLsp*fy*Dbl
                      % (usually 0.022 at ambient temp. or 0.011 at -40C)

% ============================================================================================
% ============================================================================================
% ====================              END OF INPUT DATA             ============================
% ============================================================================================
% ============================================================================================



% control parameters:

itermax    = 1000;                           % max number of iterations
ncl        = 40;                             % # of concrete layers
tolerance  = 0.001;                          % x fpc x Ag
dels       = 0.0001;                         % delta strain for default material models


% =============================================================================================




addpath('C:\CUMBIA\models')                  % directory with the user specified material models

if Ec==0
    Ec = 5000*(fpc^(1/2));                   % concrete modulus of elaticity
end

if temp < 0
    Ct = (1-0.0105*temp)*0.56*(fpc^(1/2));   % tensile strength
end

if temp >= 0
    Ct = 0.56*(fpc^(1/2));
end

eccr = Ct/Ec;                                 % concrete strain for cracking

Dsp = D-2*clb+Dh;                                     % core diameter
dcore = clb-Dh*0.5;                                   % distance to the core
P      = P*1000;                                      % axial load in Newtons
Ast    = nbl*0.25*pi*(Dbl^2);

tcl = D / ncl;                                        % thickness of concrete layers
yl  = tcl*(1:ncl);                                    % border distance conc. layer

esser = -esser;
esdam = -esdam;

switch lower(unconfined)
   case 'mu'
      [ecun,fcun] = manderun(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,'circular',D,0,0,0,0,0,dels);
   case 'mc'
      [ecun,fcun] = manderconf(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,'circular',D,0,0,0,0,0,dels,type);
   case 'mclw'
      [ecun,fcun] = manderconflw(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,'circular',D,0,0,0,0,0,dels,type);
   otherwise
    AUX = load ([unconfined,'.txt']);           % read the material model data file
    ecun = AUX(:,1)';
    fcun = AUX(:,2)';
end

switch lower(confined)
   case 'mu'
      [ec,fc] = manderun(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,'circular',D,0,0,0,0,0,dels);
   case 'mc'
      [ec,fc] = manderconf(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,'circular',D,0,0,0,0,0,dels,type);
   case 'mclw'
      [ec,fc] = manderconflw(Ec,Ast,Dh,clb,s,fpc,fyh,eco,esm,espall,'circular',D,0,0,0,0,0,dels,type); 
   otherwise
    AUX = load ([confined,'.txt']);           % read the material model data file
    ec = AUX(:,1)';
    fc = AUX(:,2)';
end


switch lower(rebar)
    case 'ks'
        [es,fs] = steelking(Es,fy,fsu,esh,esu,dels);
    case 'ra'
        [es,fs] = Raynor(Es,fy,fsu,esh,esu,dels,C1,Ey);
    otherwise
    AUX = load ([rebar,'.txt']);             % read the material model data file
    es = AUX(:,1)';
    fs = AUX(:,2)';
end
  
ecu = ec(length(ec));                         % maximum strain confined concrete
ecumander = ecu/1.5;                          % ultimate strain predicted by the original mander model

if lower(ecdam) == 'twth'
        ecdam = ecumander;
end

ec  = [-1e10 ec ec(length(ec))+dels 1e10];    % vector with strains of confined concrete
fc  = [0 fc 0 0];                             % vector with stresses of confined concrete          

ecun = [-1e10 ecun ecun(length(ecun))+dels 1e10];  % vector with strains of unconfined concrete
fcun = [0 fcun 0 0];                               % vector with stresses of unconfined concrete 

esu = es(length(es));                               % maximum strain steel
es = [es es(length(es))+dels 1e10];                % vector with strains of the steel
fs = [fs 0 0];                                     % vector with stresses of the steel
for i=1:length(es)
    esaux(i) = es(length(es)-i+1);
    fsaux(i) = fs(length(fs)-i+1);
end
es = [-esaux es(2:length(es))];                   % vector with strains of the steel
fs = [-fsaux fs(2:length(fs))];                   % vector with stresses of the steel


figure;area(ec,fc,'FaceColor','c')   
hold on
area(ecun,fcun,'FaceColor','b')
hold off;ylabel('Stress [MPa]','FontSize',16); xlabel('Strain','FontSize',16);
legend(': Confined Concrete',': Unconfined Concrete');grid on;set(gca,'Layer','top');
title('Stress-Strain Relation for Confined and Unconfined Concrete','FontSize',16)
axis([ec(2) ec(length(ec)-2) fc(1) 1.05*max(fc)]);

figure; area(es,fs,'faceColor',[0.8 0.8 0.4]); grid on
ylabel('Stress [MPa]','FontSize',16); xlabel('Strain','FontSize',16);
axis([es(3) es(length(es)-2) 1.05*fs(3) 1.05*max(fs)]);set(gca,'Layer','top');
title('Stress-Strain Relation for Reinforcing Steel','FontSize',16)

%============================== CONCRETE LAYERS ============================



yl  = sort([yl dcore D-dcore]);     % add layers to consider unconfined concrete
k = 1;

for i =1:length(yl)-1
    if yl(i)~=yl(i+1)
        yaux(k)=yl(i);
        k = k+1;
    end
end

yl = [yaux yl(length(yl))];
yc = yl-dcore;
yc = [yc(find(yc>0&yc<Dsp)), Dsp];                                           % confined concrete layers

Aaux = ((D/2)^2)*acos(1-2*yl./D)-(D/2-yl).*((D*yl-yl.^2).^(1/2));
Atc  = Aaux(1:length(Aaux))-[0,Aaux(1:length(Aaux)-1)];                     % total area of each layer

Aaux = ((Dsp/2)^2)*acos(1-2*yc./Dsp)-(Dsp/2-yc).*((Dsp*yc-yc.^2).^(1/2));
Atcc  = Aaux(1:length(Aaux))-[0,Aaux(1:length(Aaux)-1)];                    % total area of each conf. layer

k  = 1;
for i=1:length(yl)
    if yl(i)<=dcore | yl(i)>D-dcore
        conclay(i,:) = [Atc(i) 0];
    end
    if yl(i)>dcore & yl(i)<=D-dcore
        conclay(i,:) = [Atc(i)-Atcc(k) Atcc(k)];
        k = k+1;
    end
end

conclay = [yl(1)/2 0.5*(yl(1:length(yl)-1)+yl(2:length(yl)));conclay';yl]'; % [center layer|A uncon|A conf|d top layer]


%================================    REBARS     =====================================

Asb   = 0.25*pi*(Dbl^2);
r     = 0.5*( D-2*clb-Dbl);
theta = (2*pi/nbl)*[0:(nbl-1)];
distld = sort(0.5*(D-2*r)+r*sin(theta).*tan(0.5*theta));         % y coordinate of each bar
Asbs   = Asb*ones(1,nbl);                                             % area of each bar
%=============================== CORRECTED AREAS ======================================

conclay = [conclay zeros(length(conclay),1)];    % add column to store rebar area in each concrete layer 

for k=2:length(conclay)-1
    conclay(k,5)=sum(Asbs((distld<=conclay(k,4))&(distld>conclay(k-1,4))));
    if conclay(k,3)==0                                   % correct concrete areas due to rebar presence
        conclay(k,2) = conclay(k,2)-conclay(k,5);
        if conclay(k,2)<0
            disp('decrease # of layers')
            return
        end      
    else
        conclay(k,3) = conclay(k,3)-conclay(k,5);
        if conclay(k,2)<0
            disp('decrease # of layers')
            return
        end      
    end
end

%============  Define vector (def) with the deformations in the top concrete ==================

if ecu<=0.0018
    def = [0.0001:0.0001:20*ecu];  
end
if ecu>0.0018 & ecu<=0.0025
    def = [0.0001:0.0001:0.0016 0.0018:0.0002:20*ecu];
end
if ecu>0.0025 & ecu<=0.006
    def = [0.0001:0.0001:0.0016 0.0018:0.0002:0.002 0.0025:0.0005:20*ecu];
end
if ecu>0.006 & ecu<=0.012
    def = [0.0001:0.0001:0.0016 0.0018:0.0002:0.002 0.0025:0.0005:0.005 0.006:0.001:20*ecu];
end
if ecu>0.012 
    def = [0.0001:0.0001:0.0016 0.0018:0.0002:0.002 0.0025:0.0005:0.005 0.006:0.001:0.01 0.012:0.002:20*ecu];
end

np = length(def);


if P>0
    for k=1:np
        
        compch = sum(interp1(ecun,fcun,def(1)*ones(1,length(yl))).*conclay(:,2)') + ...
                 sum(interp1(ec,fc,def(1)*ones(1,length(yl))).*conclay(:,3)') + ...
                 sum(Asb*interp1(es,fs,def(1)*ones(1,length(distld))));
        if compch<P
            def = def(2:length(def));
        end
    end
end

np = length(def);

% ===============ITERATIVE PROCESS TO FIND THE MOMENT - CURVATURE RELATION: ==============================

message = 0;                                 % stop conditions

curv(1)   = 0;                               % curvatures
mom(1)    = 0;                               % moments
ejen(1)   = 0;                               % neutral axis
DF(1)     = 0;                               % force eqilibrium
vniter(1) = 0;                               % iterations
coverstrain(1) = 0;                          
corestrain(1)  = 0;                          
steelstrain(1) = 0;


tol = tolerance*0.25*pi*(D^2)*fpc;                % tolerance allowed
x(1) = D/2;
for k=1:np
    lostmomcontrol = max(mom);
    if mom(k)<(0.8*lostmomcontrol)
        message = 4; break
    end
    F = 10*tol;
    niter = 0;                                  
    while abs(F)>tol,
        niter = niter + 1;


        eec = (def(k)/x(niter))*(conclay(:,1)-(D-x(niter)));    % vector with the strains in the concrete
        ees = (def(k)/x(niter))*(distld-(D-x(niter)));          % vector with the strains in the steel

        
        fcunconf = interp1(ecun,fcun,eec);        % vector with stresses in the unconfined concr.           
        fcconf = interp1(ec,fc,eec);              % vector with stresses in the confinded concr.    
        fsteel = interp1(es,fs,ees);              % vector with stresses in the steel
        FUNCON = fcunconf.*conclay(:,2);
        FCONF  = fcconf.*conclay(:,3);
        FST    = Asb*fsteel;
        F      = sum(FUNCON) + sum(FCONF) + sum(FST) - P;
        if F>0
            x(niter+1) = x(niter) - 0.05*x(niter);
        end
        
        if F<0
            x(niter+1) = x(niter) + 0.05*x(niter);
        end
        if niter>itermax
            message = 3; break
        end
    end
    cores = (def(k)/x(niter))*abs(x(niter)-dcore);
    TF = strcmp(confined, unconfined);
    if TF == 0
        if cores>=ecu
            message = 1;break
        end
    end
    if TF == 1
        if def(k)>=ecu
            message = 1;break
        end
    end
    if abs(ees(1))>esu
        message = 2;break
    end
    ejen(1,k+1) = x(niter);    DF(1,k+1) = F;
    vniter(1,k+1) = niter;
    mom(1,k+1) = (sum(FUNCON.*conclay(:,1)) + sum(FCONF.*conclay(:,1)) + sum(FST.*distld) - P*(D/2))/(10^6);
    if mom(1,k+1)<0
        mom(1,k+1) = -0.01*mom(1,k+1);
    end
    curv(1,k+1) = 1000*def(k)/x(niter);
    coverstrain(1,k+1) = def(k);
    corestrain(1,k+1)  = cores;
    steelstrain(1,k+1) = ees(1);
    x(1) = x(niter);
    x(2:length(x))=0; 
    if message~=0
        break
    end
end


Agross = (0.25*pi*(D^2));
AsLong = nbl*Asb;
LongSteelRatio  = (nbl*Asb/(Agross));
TransvSteelRatio = pi*Dh*Dh/(s*Dsp);
AxialRatio  = (P/(fpc*Agross));

Mn     = interp1(coverstrain,mom,ecser);
esaux  = interp1(coverstrain,steelstrain,ecser);
cr = 0;  %concrete controls
    
if abs(esaux)>abs(esser) || isnan(Mn)
    Mnsteel = interp1(steelstrain,mom,esser);
    if ~isnan(Mnsteel)
        cr = 1; %steel control
        Mn = Mnsteel;
    elseif isnan(Mn) && isnan(Mnsteel)
        disp('problem with serviceability values to estimate nominal moment')
        return
    end
end

cMn = interp1(mom,ejen,Mn);

fycurvC = interp1(coverstrain, curv, 1.8*fpc/Ec);
fycurvS = interp1(steelstrain, curv, -fy/Es);

fycurv = min(fycurvC,fycurvS);                     % curvature for first yield

fyM    = interp1(curv,mom,fycurv);                  % moment for first yield

eqcurv = max((Mn/fyM)*fycurv,fycurv);

curvbilin = [0 eqcurv curv(length(curv))];
mombilin  = [0 Mn mom(length(mom))];

SectionCurvatureDuctility = curv(length(curv))/eqcurv;

figure; plot(curvbilin,mombilin,'r',curv,mom,'b--','LineWidth',2); grid on; 
xlabel('Curvature(1/m)','FontSize',16); ylabel('Moment (kN-m)','FontSize',16);     
title('Moment - Curvature Relation','FontSize',16);



for j=1:length(steelstrain)
    ffss = (-steelstrain(j)*Es);
    if ffss>fy
        ffss = fy;
    end
    Lsp(j) = kLsp*ffss*Dbl;     % Strain penetration length      
end
    


kkk = min(0.2*(fsu/fy-1),0.08);


switch lower(bending)
    case 'single'
        Lp  = max(kkk*L   + kLsp*fy*Dbl,2*kLsp*fy*Dbl);           % Plastic hinge length
        LBE = L;
    case 'double'
        Lp  = max(kkk*L/2 + kLsp*fy*Dbl,2*kLsp*fy*Dbl);            % Plastic hinge length
        LBE = L/2;
    otherwise
        disp('bending should be specified as single or double'); return;
end

% Moyer - Kowalsky Buckling model

bucritMK =0;
CuDu   = curv/eqcurv;

if SectionCurvatureDuctility>4
    
    esgr4  = -0.5*interp1(CuDu,steelstrain,4);   % resdidual growth strain at ductility 4
    escc   = 3*((s/Dbl)^(-2.5));                 % allowable steel compression strain

    for i=1:length(steelstrain)
        if CuDu(i)<1
            esgr(i) = 0;
        end
        if CuDu(i)<4 & CuDu(i)>1
            esgr(i) = (esgr4/4)*CuDu(i);
        end
        if CuDu(i)>4
            esgr(i) = -0.5*steelstrain(i);
        end
    end
    esfl = escc-esgr;

    if (-steelstrain(length(steelstrain)))>=esfl(length(esfl))
        bucritMK = 1;
        fail     = esfl-(-steelstrain);
        failCuDuMK = interp1(fail,CuDu,0);
        failesfl = interp1(fail,esfl,0);
        failss   = -interp1(fail,steelstrain,0);
        figure; plot(CuDu,-steelstrain,'-r',CuDu,esfl,'--b',failCuDuMK,failss,'mo','LineWidth',2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',12); grid on;
        xlabel('Curvature Ductility','FontSize',16); ylabel('Steel Tension Strain','FontSize',16); grid on
        legend(': Column strain ductility behavior',': Flexural Tension Strain',': Buckling')
        title('Moyer - Kowalsky Buckling Model','FontSize',16);
        
    else

        figure; plot(CuDu,-steelstrain,'-r',CuDu,esfl,'--b','LineWidth',2,'MarkerSize',6); grid on;
        xlabel('Curvature Ductility','FontSize',16); ylabel('Steel Tension Strain','FontSize',16); grid on
        legend(': Column strain ductility behavior',': Flexural Tension Strain')
        title('Moyer - Kowalsky Buckling Model','FontSize',16);
    end
end

% Berry - Eberhard Buckling model

bucritBE = 0;

if AxialRatio >= 0.30
    C0=0.006; C1=7.190; C2=3.129; C3=0.651; C4=0.227;                          % model constants
else
    C0=0.0010; C1=7.30; C2=1.30; C3=1.30; C4=3.00;                             % model constants
end

roeff = TransvSteelRatio*fyh/fpc;                                              % effective confinement ratio

rotb  = C0*(1+C1*roeff)*((1+C2*P/(Agross*fpc))^(-1))*(1+C3*LBE/D+C4*Dbl*fy/D);    % plastic rotation at the onset of bar buckling

plrot = (curv-fycurv)*(Lp)/1000;

if max(plrot)>rotb
    bucritBE = 1;
    failBE = plrot - rotb;
    failplrot  = interp1(failBE,plrot,0);
    failCuDuBE = interp1(failBE,CuDu,0);
    figure; plot(CuDu,rotb*ones(size(CuDu)),'--r',CuDu,plrot,'-b',failCuDuBE,failplrot,'mo','LineWidth',2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',12); grid on;
    xlabel('Curvature Ductility','FontSize',16); ylabel('Plastic Rotation','FontSize',16); grid on
    legend(': Plastic Rotation for Buckling',': Plastic Rotation',': Buckling')
    title('Berry - Eberhard Buckling Model','FontSize',16);
else

    figure; plot(CuDu,rotb*ones(size(CuDu)),'--r',CuDu,plrot,'-b','LineWidth',2, 'MarkerEdgeColor','k',...
        'MarkerFaceColor','g','MarkerSize',12); grid on;
    xlabel('Curvature Ductility','FontSize',16); ylabel('Plastic Rotation','FontSize',16); grid on
    legend(': Plastic Rotation for Buckling',': Plastic Rotation')
    title('Berry - Eberhard Buckling Model','FontSize',16);
end



% Flexure deflection:

switch lower(bending)
    case 'single'
        for i=1:length(curv)
            if coverstrain(i)<eccr
                displf(i) = curv(i)*((L/1000)^2)/3;
            end
            if coverstrain(i)>eccr & curv(i)<fycurv
                displf(i) = curv(i) * (((L+Lsp(i))/1000)^2)/3;
            end
            if curv(i)>=fycurv
                displf(i) = (curv(i)-fycurv*(mom(i)/fyM))*(Lp/1000)*((L+Lsp(i)-0.5*Lp)/1000) +...
                            (fycurv*(((L+Lsp(i))/1000)^2)/3)*(mom(i)/fyM);
            end
        end
        Force = mom/(L/1000);
    case 'double'
        for i=1:length(curv)
            if coverstrain(i)<eccr
                displf(i) = curv(i)*((L/1000)^2)/6;
            end
            if coverstrain(i)>eccr & curv(i)<fycurv
                displf(i) = curv(i) * (((L+2*Lsp(i))/1000)^2)/6;
            end
            if curv(i)>=fycurv
                displf(i) = (curv(i)-fycurv*(mom(i)/fyM))*(Lp/1000)*((L+2*(Lsp(i)-0.5*Lp))/1000) +...
                            (fycurv*(((L+2*Lsp(i))/1000)^2)/6)*(mom(i)/fyM);
            end
        end
        Force = 2*mom/(L/1000);
    otherwise
        disp('bending should be specified as single or double'); return;
end


% Shear deflection:

G     = 0.43*Ec;
As    = 0.9*Agross;
Ig    = pi*(D^4)/64;
Ieff  = (Mn*1000/(Ec*(10^6)*eqcurv))*(10^12);

beta  = min(0.5+20*LongSteelRatio,1);

switch lower(bending)
    case 'single'
        alpha = min(max(1,3-L/D),1.5);
    case 'double'
        alpha = min(max(1,3-L/(2*D)),1.5);
end

Vc1   = 0.29*alpha*beta*0.8*(fpc^(1/2))*Agross/1000;

kscr  = ((0.39*TransvSteelRatio)*0.25*Es*((0.8*D/1000)^2)/(0.25+10*(0.39*TransvSteelRatio)))*1000;

switch lower(bending)
    case 'single'
        ksg   = (G*As/L)/1000;
        kscr  = (kscr/L);
        forcebilin = mombilin/(L/1000);
    case 'double'
        ksg   = (G*As/(L/2))/1000;
        kscr  = (kscr/(L/2));
        forcebilin = 2*mombilin/(L/1000);
end
kseff = ksg*(Ieff/Ig);
aux = (Vc1/kseff)/1000;
aux2 = 0;
momaux = mom;
for i=1:length(curv)
    if momaux(i)<=Mn & Force(i)<Vc1
        displsh(i) = (Force(i)/kseff)/1000;
    end
    if momaux(i)<=Mn & Force(i)>=Vc1
        displsh(i) = ((Force(i)-Vc1)/kscr)/1000+aux;
    end
    if momaux(i)>Mn
        momaux = 4*momaux;
        aux3=i-aux2;
        aux2 = aux2 + 1;
        displsh(i) = (displf(i)/displf(i-1))*displsh(i-1);
    end
end

displ = displsh + displf;

% bilinear approx:

dy1 = interp1(curv,displ,fycurv);
dy  = (Mn/fyM)*dy1;
du  = displ(length(displ));
displbilin  = [0 dy du];
Dduct = displ/dy;
DisplDuct = max(Dduct);

dy1f = interp1(curv,displf,fycurv);
dyf  = (Mn/fyM)*dy1f;


% Shear Strength:

Vs    = (0.5*pi*(0.25*pi*(Dh^2))*fyh*cot(pi/6)*(D-clb+0.5*Dh-cMn)/s)/1000;
Vsd   = (0.5*pi*(0.25*pi*(Dh^2))*fyh*cot((35/180)*pi)*(D-clb+0.5*Dh-cMn)/s)/1000;
beta  = min(0.5+20*LongSteelRatio,1);
Dductf = displ/dyf;

switch lower(bending)
    case 'single'
        alpha = min(max(1,3-L/D),1.5);
        if P>0
            Vp = (P*(D-cMn)/(2*L))/1000;
        else
            Vp = 0;
        end

    case 'double'
        alpha = min(max(1,3-L/(2*D)),1.5);
        if P>0
            Vp = (P*(D-cMn)/(L))/1000;
        else
            Vp = 0;
        end

end

switch lower(ductilitymode)
    case 'uniaxial'
        for i =1:length(Dductf)
            Vc(i) = alpha*beta*min(max(0.05,0.37-0.04*Dductf(i)),0.29)*0.8*(fpc^(1/2))*Agross/1000;
        end
    case 'biaxial'
        for i =1:length(Dductf)
            Vc(i) = alpha*beta*min(max(0.05,0.33-0.04*Dductf(i)),0.29)*0.8*(fpc^(1/2))*Agross/1000;
        end
end

Vcd = 0.862*Vc;
Vpd = 0.85*Vp;
V   = Vc + Vs + Vp;
Vd  = 0.85 * ( Vcd +Vsd + Vpd);
criteria = 1;
if V(length(V))<Force(length(Force))
    failure   = V-Force;
    faildispl = interp1(failure,displ,0);
    failforce = interp1(displ,Force,faildispl);
    failduct  = interp1(displ,Dduct,faildispl);
    failmom   = interp1(displ,mom,faildispl);
    failcurv  = interp1(displ,curv,faildispl);
    failCuDu  = interp1(displ,CuDu,faildispl);
    switch lower(bending)
        case 'single'
            if faildispl<=2*dy
                criteria = 2;
            elseif faildispl<8*dy
                criteria = 3;
            else
                criteria = 4;
            end
        case 'double'
            if faildispl<=1*dy
                criteria = 2;
            elseif faildispl<7*dy
                criteria = 3;
            else
                criteria = 4;
            end
    end
end



Ieq = (Mn/(eqcurv*Ec))/1000; % equivalent I for NL THA
Bi = 1/(((mombilin(2))/(curvbilin(2)))/((mombilin(3)-mombilin(2))/(curvbilin(3)-curvbilin(2)))); % Bilinear factor


% Limit States:

displdam = 0; displser = 0; Dductdam = 0; Dductser = 0;
curvdam = 0;   curvser  = 0; CuDudam  = 0; CuDuser = 0;
coverstraindam = 0; coverstrainser = 0;
steelstraindam = 0; steelstrainser = 0;
momdam = 0; momser = 0; Forcedam = 0; Forceser = 0; 

if max(coverstrain) > ecser | max(abs(steelstrain)) > abs(esser)
    
    if max(coverstrain) > ecdam | max(abs(steelstrain)) > abs(esdam)

        displdamc = interp1(coverstrain,displ,ecdam);
        displdams = interp1(steelstrain,displ,esdam);
        displdam  = min (displdamc,displdams);
        Dductdam  = interp1(displ,Dduct,displdam);
        curvdam   = interp1(displ,curv,displdam);
        CuDudam   = interp1(displ,CuDu,displdam);
        coverstraindam = interp1(displ,coverstrain,displdam);
        steelstraindam = interp1(displ,steelstrain,displdam);
        momdam   = interp1(displ,mom,displdam);
        Forcedam = interp1(displ,Force,displdam);
        
    end

        displserc = interp1(coverstrain,displ,ecser);
        displsers = interp1(steelstrain,displ,esser);
        displser  = min (displserc,displsers);
        Dductser  = interp1(displ,Dduct,displser);
        curvser   = interp1(displ,curv,displser);
        CuDuser   = interp1(displ,CuDu,displser);
        coverstrainser = interp1(displ,coverstrain,displser);
        steelstrainser = interp1(displ,steelstrain,displser);
        momser   = interp1(displ,mom,displser);
        Forceser = interp1(displ,Force,displser);

end

outputlimit = [ coverstrainser steelstrainser momser Forceser curvser CuDuser displser Dductser
                coverstraindam steelstraindam momdam Forcedam curvdam CuDudam displdam Dductdam
                max(coverstrain) min(steelstrain) mom(length(mom)) Force(length(Force)) max(curv) max(CuDu) max(displ) max(Dduct) ];


figure; 
plot(displbilin,forcebilin,'b--','LineWidth',2, 'DisplayName',': bilinear approximation'); grid on;
xlabel('Displacement(m)','FontSize',16); ylabel('Force (kN)','FontSize',16);
title('Force - Displacement Relation','FontSize',16);
hold on

plot(displ,Force,'k','LineWidth',2, 'DisplayName',': total response');
plot(displ,V,'r:','LineWidth',2, 'DisplayName',': shear capacity (assessment)');
plot(displ,Vd,'m:','LineWidth',2, 'DisplayName',': shear capacity (design)');


if criteria ~=1
    plot(faildispl,failforce,'mo','LineWidth',2, 'DisplayName',': shear failure','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
end
if bucritMK == 1
    buckldispl = interp1(CuDu,displ,failCuDuMK);
    bucklforce = interp1(CuDu,Force,failCuDuMK);
    plot(buckldispl,bucklforce,'k*','LineWidth',2, 'DisplayName',': buckling (M & K)','MarkerEdgeColor','k','MarkerSize',10,'MarkerFaceColor','g')
end
if bucritBE==1
    buckldisplBE = interp1(CuDu,displ,failCuDuBE);
    bucklforceBE = interp1(CuDu,Force,failCuDuBE);
    plot(buckldisplBE,bucklforceBE,'ks','LineWidth',2, 'DisplayName',': buckling (B & E)','MarkerEdgeColor','k','MarkerSize',10,'MarkerFaceColor','g')
end

legend('Location','southeast')

hold off

%==========================================================================

pointsdam = find(displ<=displdam);
pointsser = find(displ<=displser);

figure;
area(displ,Force,'FaceColor',[0.4 0.8 0.6], 'DisplayName',': ultimate zone');grid on; hold on;
area(displ(pointsdam),Force(pointsdam),'FaceColor',[0.4 0.6 0.6], 'DisplayName',': damage control zone');
area(displ(pointsser),Force(pointsser),'FaceColor',[0.4 0.4 0.6], 'DisplayName',': serviceability zone');

plot(displbilin,forcebilin,'b--','LineWidth',2, 'DisplayName',': bilinear approximation');set(gca,'Layer','top');
plot(displ,Force,'k','LineWidth',2, 'DisplayName',': total response');set(gca,'Layer','top');set(gca,'Layer','top');
plot(displ,V,'r:','LineWidth',2, 'DisplayName',': shear capacity (assessment)');set(gca,'Layer','top');
plot(displ,Vd,'m:','LineWidth',2, 'DisplayName',': shear capacity (design)');set(gca,'Layer','top');

xlabel('Displacement(m)','FontSize',16); ylabel('Force (kN)','FontSize',16);
title('Potential Deformation Limit States','FontSize',16);

if criteria ~=1
    plot(faildispl,failforce,'mo','LineWidth',2, 'DisplayName',': shear failure','MarkerEdgeColor','k','MarkerFaceColor','g','MarkerSize',10)
end
if bucritMK == 1
    buckldispl = interp1(CuDu,displ,failCuDuMK);
    bucklforce = interp1(CuDu,Force,failCuDuMK);
    plot(buckldispl,bucklforce,'k*','LineWidth',2, 'DisplayName',': buckling (M & K)','MarkerEdgeColor','k','MarkerSize',10,'MarkerFaceColor','g')
end
if bucritBE==1
    buckldisplBE = interp1(CuDu,displ,failCuDuBE);
    bucklforceBE = interp1(CuDu,Force,failCuDuBE);
    plot(buckldisplBE,bucklforceBE,'ks','LineWidth',2, 'DisplayName',': buckling (B & E)','MarkerEdgeColor','k','MarkerSize',10,'MarkerFaceColor','g')
end

legend('Location','northwest')

hold off

%==========================================================================

output = [coverstrain;corestrain;ejen;steelstrain;mom;curv;Force;displsh;displf;displ;V;Vd];
outputbilin = [curvbilin;mombilin;displbilin;forcebilin];

Acore = 0.25*pi*(Dsp^2);

PCid   = interp1(ec,fc,csid)*(Acore-AsLong)+interp1(ecun,fcun,csid)*(Agross-Acore)+AsLong*interp1(es,fs,csid); % compression force for yield surface
PTid   = AsLong*interp1(es,fs,ssid);                                                                           % tensile force for yield surface


fid = fopen([name,'.xls'],'w');
fprintf(fid, '\n\nCircular Section\n\n');
switch lower(confined)
    case 'mclw'
      fprintf(fid,'lightweight concrete\n');
    otherwise
      fprintf(fid,'normalweight concrete\n');
end
fprintf(fid, 'Diameter:  %5.1f mm\n',D);
fprintf(fid, 'cover to longitudinal bars:  %4.1f mm\n',clb);
fprintf(fid, 'number of longitudinal bars:  %4.0f\n',nbl);
fprintf(fid, 'diameter of longitudinal bars:  %4.1f mm\n',Dbl);
fprintf(fid, 'diameter of transverse steel:  %4.1f mm\n',Dh);
fprintf(fid, 'spacing of transverse steel:  %4.1f mm\n',s);
switch lower(type)
    case 'spirals'
      fprintf(fid,'type of tranverse reinforcement: spirals\n');
    case 'hoops'
      fprintf(fid,'type of tranverse reinforcement: hoops\n');
end
fprintf(fid, 'axial load:  %8.2f kN\n',P/1000);
fprintf(fid, 'concrete compressive strength:  %3.2f MPa\n',fpc);
fprintf(fid, 'long steel yielding stress:  %4.2f MPa\n',fy);
fprintf(fid, 'long steel max. stress:  %4.2f MPa\n',max(fs));
fprintf(fid, 'transverse steel yielding stress:  %4.2f MPa\n',fyh);
fprintf(fid, 'Member Length:  %5.1f mm\n',L);
switch lower(bending)
    case 'single'
      fprintf(fid,'Single Bending\n');
    case 'double'
      fprintf(fid,'Double Bending\n');
end
switch lower(ductilitymode)
    case 'uniaxial'
      fprintf(fid,'Uniaxial Bending\n');
    case 'biaxial'
      fprintf(fid,'Biaxial Bending\n');
end
fprintf(fid, 'Longitudinal Steel Ratio:  %1.3f\n',LongSteelRatio);
fprintf(fid, 'Transverse Steel Ratio:  %1.3f\n',TransvSteelRatio);
fprintf(fid, 'Axial Load Ratio:  %1.3f\n\n',AxialRatio);
fprintf(fid, 'Cover\tCore\tN.A\tSteel\tMoment\tCurvature\tForce\tSh displ.\tFl displ.\tTotal displ.\tShear(assess.)\tShear(design)\n'); 
fprintf(fid, 'Strain\tStrain\t[mm]\tStrain\t[kN-m]\t[1/m]\t[kN]\t[m]\t[m]\t[m]\t[kN]\t[kN]\n'); 
fprintf(fid, '%1.5f\t%1.5f\t%4.2f\t%1.5f\t%8.2f\t%1.5f\t%8.2f\t%1.5f\t%1.5f\t%1.5f\t%8.2f\t%8.2f\n',output);
fprintf(fid, ' \n');
fprintf(fid, 'Bilinear Approximation:\n\n');
fprintf(fid, 'Curvature\tMoment\tDispl.\tForce\n'); 
fprintf(fid, '[1/m]\t[kN-m]\t[m]\t[kN]\n'); 
fprintf(fid,'%1.5f\t%8.2f\t%1.5f\t%8.2f\n',outputbilin);
fprintf(fid, ' \n');
switch message
   case 1
      fprintf(fid,' *** concrete strain exceeds maximum ***');
   case 2
      fprintf(fid,' *** steel strain exceeds maximum ***');
   case 3
      fprintf(fid,' *** number of iteration exceeds maximum ***');
   case 4
      fprintf(fid,' *** excessive lost of strength ***');
end
fprintf(fid, '\n\nMoment for First Yielding:  %8.2f kN-m\n',fyM);
fprintf(fid, 'Curvature for First Yielding:  %1.5f 1/m\n',fycurv);
fprintf(fid, 'Potential Section Nominal Moment:  %8.2f kN-m\n',Mn);
fprintf(fid, 'Equivalent Curvature:  %1.5f 1/m\n',eqcurv);
fprintf(fid, 'Potential Section Curvature Ductility:  %3.2f\n',SectionCurvatureDuctility);
fprintf(fid, 'Potential Displacement Ductility:  %3.2f\n',DisplDuct);
fprintf(fid, ' \n');
switch criteria
    case 1
        fprintf(fid,' *** flexural failure ***');
    case 2
        fprintf(fid,' *** brittle shear failure ***');
        fprintf(fid, '\n\nDisplacement for Shear Failure:  %1.5f m\n',faildispl);
        fprintf(fid, 'Displacement Ductility at Shear Failure:  %8.2f\n',failduct);
        fprintf(fid, 'Force for Shear Failure:  %8.2f kN\n',failforce);
        fprintf(fid, 'Curvature for Shear Failure:  %1.5f 1/m\n',failcurv);
        fprintf(fid, 'Curvature Ductility at Shear Failure:  %8.2f\n',failCuDu);
        fprintf(fid, 'Moment for Shear Failure:  %8.2f kN-m\n',failmom);
        
    case 3
        fprintf(fid,' *** shear failure at some ductility ***');
        fprintf(fid, '\n\nDisplacement for Shear Failure:  %1.5f m\n',faildispl);
        fprintf(fid, 'Displacement Ductility at Shear Failure:  %8.2f\n',failduct);
        fprintf(fid, 'Force for Shear Failure:  %8.2f kN\n',failforce);
        fprintf(fid, 'Curvature for Shear Failure:  %1.5f 1/m\n',failcurv);
        fprintf(fid, 'Curvature Ductility at Shear Failure:  %8.2f\n',failCuDu);
        fprintf(fid, 'Moment for Shear Failure:  %8.2f kN-m\n',failmom);
        
    case 4
        fprintf(fid,' *** ductil shear failure ***');
        fprintf(fid, '\n\nDisplacement for Shear Failure:  %1.5f m\n',faildispl);
        fprintf(fid, 'Displacement Ductility at Shear Failure:  %8.2f\n',failduct);
        fprintf(fid, 'Force for Shear Failure:  %8.2f kN\n',failforce);
        fprintf(fid, 'Curvature for Shear Failure:  %1.5f 1/m\n',failcurv);
        fprintf(fid, 'Curvature Ductility at Shear Failure:  %8.2f\n',failCuDu);
        fprintf(fid, 'Moment for Shear Failure:  %8.2f kN-m\n',failmom);
        
end
fprintf(fid, ' \n');
if bucritMK == 1
    bucklDd =   interp1(CuDu,Dduct,failCuDuMK);
    bucklcurv = interp1(CuDu,curv,failCuDuMK);
    bucklmom  = interp1(CuDu,mom,failCuDuMK);
    fprintf(fid,'Moyer - Kowalsky buckling model:\n');
    fprintf(fid, '\nCurvature Ductility for Buckling:  %8.2f\n',failCuDuMK);
    fprintf(fid, 'Curvature at Buckling:  %3.5f m\n',bucklcurv);
    fprintf(fid, 'Displacement Ductility at Buckling:  %8.2f\n',bucklDd);
    fprintf(fid, 'Displacement at Buckling:  %3.5f m\n',buckldispl);
    fprintf(fid, 'Force for Buckling:  %8.2f kN\n',bucklforce);
    fprintf(fid, 'Moment for Buckling:  %8.2f kN\n',bucklmom);
    
end
fprintf(fid, '\n');
if bucritBE == 1
    bucklDdBE = interp1(CuDu,Dduct,failCuDuBE);
    bucklcurvBE = interp1(CuDu,curv,failCuDuBE);
    bucklmomBE  = interp1(CuDu,mom,failCuDuBE);
    fprintf(fid,'Berry - Eberhard buckling model:\n');
    fprintf(fid, '\nCurvature Ductility for Buckling:  %8.2f\n',failCuDuBE);
    fprintf(fid, 'Curvature at Buckling:  %3.5f m\n',bucklcurvBE) 
    fprintf(fid, 'Displacement Ductility at Buckling:  %8.2f\n',bucklDdBE);
    fprintf(fid, 'Displacement at Buckling:  %3.5f m\n',buckldisplBE);
    fprintf(fid, 'Force for Buckling:  %8.2f kN\n',bucklforceBE);
    fprintf(fid, 'Moment for Buckling:  %8.2f kN\n',bucklmomBE);
    
end
fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '== Potential Deformation Limit States (serviceability/damage control/ultimate) ==\n\n');
fprintf(fid, 'Cover\tSteel\tMoment\tForce\tCurvature\tCurvature\tDisplacement\tDisplacement\n'); 
fprintf(fid, 'Strain\tStrain\t[kN-m]\t[kN]\t[1/m]\tDuctility\t[m]\tDuctility\n'); 
fprintf(fid, '%1.5f\t%1.5f\t%8.2f\t%8.2f\t%1.5f\t%3.2f\t%2.5f\t%3.2f\n',outputlimit');
fprintf(fid, '\nDeformation Limit States Citeria :\n');
fprintf(fid, '   serviceability concrete strain: %1.4f\n',ecser); 
fprintf(fid, '   serviceability steel strain: %1.4f\n',esser);
fprintf(fid, '   damage control concrete strain: %1.4f\n',ecdam);
fprintf(fid, '   damage control steel strain: %1.4f\n',esdam);

switch lower(confined)
   case 'mc'
      fprintf(fid,'\nOriginal Mander Model Ultimate Concrete Strain: %1.4f\n',ecumander);
end


fprintf(fid, '\nfor non-linear THA:\n\n');
fprintf(fid, 'E: %10.2f Pa\n',Ec*(10^6));
fprintf(fid, 'G: %10.2f Pa\n', G*(10^6));
fprintf(fid, 'A: %10.4f m2\n', Agross/(10^6));
fprintf(fid, 'I: %10.6f m4\n', Ieq);
fprintf(fid, 'Bi-Factor: %1.3f\n', Bi);
fprintf(fid, 'Hinge Length: %1.3f m\n', Lp/1000);
fprintf(fid, 'Tension Yield: %10.2f N\n', PTid);
fprintf(fid, 'Compression Yield: %10.2f N\n', PCid);
fprintf(fid, 'Moment Yield: %10.2f N-m\n', Mn*1000);


if interaction ~= 'y'
    fclose(fid);
    return
end





PP  = [-0.90*PTid:0.30*PTid:0 0.05*fpc*Agross:0.05*fpc*Agross:0.7*PCid];               % vector with axial loads for interaction diagram
nPP = length(PP);

for i=1:nPP
    
    curv=0;mom=0;ejen=0;DF=0;vniter=0;coverstrain=0;corestrain=0;
    steelstrain=0;x=0;message=0;
    
    if ecu<=0.0018
        def = [0.0001:0.0001:20*ecu];
    end
    if ecu>0.0018 & ecu<=0.0025
        def = [0.0001:0.0001:0.0016 0.0018:0.0002:20*ecu];
    end
    if ecu>0.0025 & ecu<=0.006
        def = [0.0001:0.0001:0.0016 0.0018:0.0002:0.002 0.0025:0.0005:20*ecu];
    end
    if ecu>0.006 & ecu<=0.012
        def = [0.0001:0.0001:0.0016 0.0018:0.0002:0.002 0.0025:0.0005:0.005 0.006:0.001:20*ecu];
    end
    if ecu>0.012
        def = [0.0001:0.0001:0.0016 0.0018:0.0002:0.002 0.0025:0.0005:0.005 0.006:0.001:0.01 0.012:0.002:20*ecu];
    end

    np = length(def);

    if PP(i)>0
        for j=1:np
                compch = sum(interp1(ecun,fcun,def(1)*ones(1,length(yl))).*conclay(:,2)') + ...
                sum(interp1(ec,fc,def(1)*ones(1,length(yl))).*conclay(:,3)') + ...
                sum(Asb*interp1(es,fs,def(1)*ones(1,length(distld))));
            if compch<PP(i)
                def = def(2:length(def));
            end
        end
    end
    
    np = length(def);

    curv(1)   = 0;                               % curvatures
    mom(1)    = 0;                               % moments
    ejen(1)   = 0;                               % neutral axis
    DF(1)     = 0;                               % force eqilibrium
    vniter(1) = 0;                               % iterations
    coverstrain(1) = 0;
    corestrain(1)  = 0;
    steelstrain(1) = 0;
    tol = tolerance*0.25*pi*(D^2)*fpc;;                % tolerance allowed
    x(1) = D/2;
    
    for k=1:np
        F = 10*tol;
        niter = 0;
        while abs(F)>tol,
            niter = niter + 1;

            eec = (def(k)/x(niter))*(conclay(:,1)-(D-x(niter)));    % vector with the strains in the concrete
            ees = (def(k)/x(niter))*(distld-(D-x(niter)));          % vector with the strains in the steel

            fcunconf = interp1(ecun,fcun,eec);        % vector with stresses in the unconfined concr.
            fcconf = interp1(ec,fc,eec);              % vector with stresses in the confinded concr.
            fsteel = interp1(es,fs,ees);              % vector with stresses in the steel
            FUNCON = fcunconf.*conclay(:,2);
            FCONF  = fcconf.*conclay(:,3);
            FST    = Asb*fsteel;
            F      = sum(FUNCON) + sum(FCONF) + sum(FST) - PP(i);
            if F>0
                x(niter+1) = x(niter) - 0.05*x(niter);
            end

            if F<0
                x(niter+1) = x(niter) + 0.05*x(niter);
            end
            if niter>itermax
                message = 3; break
            end
        end
        cores = (def(k)/x(niter))*abs(x(niter)-dcore);
        TF = strcmp(confined, unconfined);
        if TF == 0
            if cores>=ecu
                message = 1;break
            end
        end
        if TF == 1
            if def(k)>=ecu
                message = 1;break
            end
        end
        if abs(ees(1))>esu
            message = 2;break
        end
        ejen(1,k+1) = x(niter);    DF(1,k+1) = F;
        vniter(1,k+1) = niter;
        mom(1,k+1) = (sum(FUNCON.*conclay(:,1)) + sum(FCONF.*conclay(:,1)) + sum(FST.*distld) - PP(i)*(D/2))/(10^6);
        curv(1,k+1) = 1000*def(k)/x(niter);
        coverstrain(1,k+1) = def(k);
        corestrain(1,k+1)  = cores;
        steelstrain(1,k+1) = ees(1);
        x(1) = x(niter);
        x(2:length(x))=0;
        if message~=0
            break
        end
    end
     Mni(i)     = interp1(coverstrain,mom,csid);
    esaux      = interp1(coverstrain,steelstrain,csid);
    cr = 0;  %concrete control
    
    if abs(esaux)>abs(ssid) || isnan(Mni(i))
        Mnsteel     = interp1(steelstrain,mom,-ssid);
        if ~isnan(Mnsteel)
            cr = 1; %steel control
            Mni(i) = Mnsteel;
        elseif isnan(Mni(i)) && isnan(Mnsteel)
            disp('problem with strain values to estimate PM interaction at axial force:')
            disp(PP(i))
            return
        end
        
    end
    if cr==0
        eci(i)  = csid;
        esi(i)  = esaux;
    end
    if cr==1
        esi(i) = -ssid;
        eci(i) = interp1(steelstrain,coverstrain,-ssid);
    end
    mess(i)  = message;
end

Mni = [0 Mni 0];

PPn  = [-PTid  PP PCid];


MB = max(Mni);
PB = PPn(find(Mni==MB));

PB13 = (1/3)*PB;
MB13 = interp1(PPn,Mni,PB13);

PB23 = (2/3)*PB;
MB23 = interp1(PPn,Mni,PB23);

MB0 = interp1(PPn,Mni,0);

PPL = [-PTid 0 PB13 PB23 PB PCid];
MnL = [0 MB0 MB13 MB23 MB 0];

outputi  = [Mni;PPn/1000];

figure; plot(Mni,PPn/1000,'-ro',MnL,PPL/1000,'--bs','LineWidth',2,'MarkerSize',6); grid on;
        xlabel('Moment [kN-m]','FontSize',16); ylabel('Axial Load [kN]','FontSize',16);
        legend(': Interaction Diagram',': Approximation for NLTHA');
        title('Interaction Diagram','FontSize',16);
        
       
        
fprintf(fid, '\n\n=*=*=*=*=*=*=*=*=**=*=*=*=*=*=**=*=*\n\n');
fprintf(fid, '\n\nInteraction Surface\n\n');
fprintf(fid, 'Concrete limit strain:  %1.4f\n',csid);
fprintf(fid, 'Steel limit strain:  %1.4f\n',ssid);
fprintf(fid, '\nMoment\tAxial Load\n'); 
fprintf(fid, '[kN-m]\t[kN]\n'); 
fprintf(fid, '%8.2f\t%8.2f\n',outputi);
fprintf(fid, ' \n');
fprintf(fid, 'NLTHA Approximation:\n\n');
fprintf(fid, 'PT:  %6.1f kN\n',-PTid/1000);
fprintf(fid, 'PC:  %6.1f kN\n',PCid/1000);
fprintf(fid, '     PB:  %6.1f kN        MB:  %6.1f kN-m\n',PB/1000,MB);
fprintf(fid, '(1/3)PB:  %6.1f kN    (1/3)MB:  %6.1f kN-m\n',PB13/1000,MB13);
fprintf(fid, '(2/3)PB:  %6.1f kN   (2/3)MB:  %6.1f kN-m\n',PB23/1000,MB23);

fclose(fid);



figure; plot(PPn(2:i+1)/1000,eci,'--bs','LineWidth',2,'MarkerSize',6); grid on;
        xlabel('Axial Force [kN]','FontSize',16); ylabel('Concrete Strain','FontSize',16);
        
figure; plot(PPn(2:i+1)/1000,esi,'--bs','LineWidth',2,'MarkerSize',6); grid on;
        xlabel('Axial Force [kN]','FontSize',16); ylabel('Steel Strain','FontSize',16);
  
    

  