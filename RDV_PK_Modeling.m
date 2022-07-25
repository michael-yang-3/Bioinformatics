% Michael Yang | my2699
% Due Oct 04, 2020
% BMEN E4530
%% Homework 2 Problem 5

% to run your code simply hit "Run" or use the shortcut "command(Mac)/
% control(Windows)+enter"

clear; clc % this is just good practice, to clear your variables and 
% command window each time you run your code

% BEGIN HERE
% set time span, in days
tspan = 1;


% volumes, in L
V_1 = 2.8*10^(-3);
V_2 = 1.1*10^(-7);
V_3 = 1.7*10^(-5);
% V_4 = 3.5*10^(-3); Not necessary 
V_5 = 5.8*10^(-5); % Using V4T as it is necessary for NTP


% initial conditions
A1_0 = 1.12*10^(-4); %grams, multiplied init_conc by V1 converted to mL
% Rest should be 0, as stated on Piazza at time 0
A2_0 = 0; 
A3_0 = 0; 
A4_0 = 0; 
A5_0 = 0; 
A6_0 = 0; 
A7_0 = 0; 
y0 = [A1_0 A2_0 A3_0 A4_0 A5_0 A6_0 A7_0];


% molecular weights, in g/mol
RDV = 602.6;
Alanine_metabolite = 442.3;
Parent_nucleoside = 291.3;
NTP = 527.2;


% set relative and absolute tolerances for increased accuracy
opts = odeset('RelTol', 1*10^(-9), 'AbsTol', 1*10^(-10));


% run ode45
% Where A1T = A(4), A2T = A(5), A3T = A(6)...
[time, S] = ode15s(@odefun, [0 tspan], y0, opts);
%ode45 created a really weird graph 2, so I used ode15s instead


% Convert into concentration plots for S 1-4
% Divide what we get in A by MW (to get mols) then divide by volume for
% concentration and convert mol/L into umol/L
S(:,1) = ((S(:,1) / RDV) / V_1)*(10^6); % plasma RDV
S(:,2) = ((S(:,2) / Alanine_metabolite) / V_2)*(10^6); % plasma Alanine Metabolite
S(:,3) = ((S(:,3) / Parent_nucleoside) / V_3)*(10^6); % plasma Parent Nucleoside
S(:,7) = ((S(:,7) / NTP) / V_5)*(10^6);% compartment for NTP in PBMC


% plot concentration vs. time graphs (uM vs. days)
% plot(time, S(:,1));
% please use subplot() to get all graphs in same figure
figure;
subplot(1, 4, 1);
plot(time, S(:,1));
title('Plasma RDV Concentration Over Time');
xlabel('Time(days)');
ylabel('Plasma GS5734 (uM)');
subplot(1, 4, 2);
plot(time, S(:,2));
title('Plasma Alanine-metabolite Concentration');
xlabel('Time(days)');
ylabel('Plasma ALA-MET (uM)');
subplot(1, 4, 3);
plot(time, S(:,3));
title('Plasma Parent Nucleoside Concentration');
xlabel('Time(days)');
ylabel('Plasma GS-441524 (uM)');
subplot(1, 4, 4);
plot(time, S(:,7));
title('Concentration of NTP in PBMC');
xlabel('Time(days)');
ylabel('NTP in PBMC (uM)');

% create function for system of ODEs titled "odefun"
%%
function dAdt = odefun(~,A)
% rate constants (assume it is in units of days^-1)
    k3e = 1.0; % best guess value
    kc3 = log(2)/(1/24); % best guess value
    kc1 = log(2)/(1/24);
    kc2 = log(2)/1;
    kcT4 = log(2)/1;
    k1T = 33.0;
    k2T = 5.86*10^3;
    k1e = 0.02;
    k2e = 12.3;
    k12 = 1.0;
    k23 = 989.6;
    k34 = 158.4;
    k3T = 7.9;
    
    
% create 7 x 1 matrix of zeros for initialization
% init = [0; 0; 0; 0; 0; 0; 0];

    
% system of differential equations (use A as your main variable)
% Where A1T = A(4), A2T = A(5), A3T = A(6)...
dAdt(1) = -kc1*A(1) - k12*A(1) - k1T*A(1) + k1e*A(4); %A1
dAdt(2) = -kc2*A(2) + k12*A(1) - k23*A(2) - k2T*A(2) + k2e*A(5); %A2
dAdt(3) = -kc3*A(3) + k23*A(2) - k3T*A(3) + k3e*A(6); %A3
dAdt(4) = k1T*A(1) - k12*A(4) - k1e*A(4); %A1T
dAdt(5) = k2T*A(2) + k12*A(4) - k23*A(5) - k2e*A(5); %A2T
dAdt(6) = k3T*A(3) + k23*A(5) - k34*A(6) - k3e*A(6); %A3T
dAdt(7) = k34*A(6) - kcT4*A(7); %A4T

dAdt = dAdt(:); %createa  column vector, found on mathworks forums

end