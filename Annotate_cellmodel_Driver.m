%% 
% This file initialises parameters and calls the solver routine.
%
% NOTE 1: 'placeholder' variables can be used to simulate effect of a
% heterologous protein.
%
% NOTE 2: z (zombie) variables have been removed for clarity. These are
% only required for parameterisation, and can be manually added back.
%
% NOTE 3: For extra development, see the arXiv paper from Diego Oyarzun
% on how to apply this model to other SynBio designs:
% "Prediction of cellular burden with host-circuit models"
% https://arxiv.org/abs/2004.00995



%% USEFUL PLOTTING COMPONENTS

%%% Useful global settings
set(0, 'DefaultAxesFontSize', 22) % Axis font size
set(0,'DefaultLineLineWidth', 2); % For line plots
set(0,'DefaultLineMarkerSize', 10) % For scatter plots

%%% Useful colour defaults
LightGrey = [0.9 0.9 0.9];
MediumGrey = [0.55 0.55 0.55];
DarkGrey = [0.2 0.2 0.2];
MATLABblue = [3 115 189]/255;
MATLABorange = [217 87 30]/255;

%%% Additional commands when plotting (add after plot function)
% axis square
% grid on



%% PARAMETERS AND CONSTANTS

% Parameters - cut to 4 sig fig where necessary
% * - value inferred through parameter fitting

% Define parameters
s0 = 10000;         % External nutrient                 [molecs]
ns = 0.5;           % Nutrient efficiency               []
nr = 7549;          % Ribosome length                   []
nx = 300;           % Length of non-ribosomal proteins  [aa/molecs]
gmax = 1260;        % Max. trans. elongation rate       [aa/min/molecs]
Kp = 180.1;         % = gmax/Kgamma = 1260/7            [aa cell/min/molecs^2] *
vet = 726;          % Max. nutrient import threshold    [/min]
Ket = 1000;         % Nutrient import threshold         [molecs]
vem = 5800;         % Max. enzymatic rate               [/min]
Kem = 1000;         % Enzymatic threshold               [molecs/cell]
wr = 930.0;         % Max. ribosome transcr. rate       [molecs/min/cell] *
we = 4.139;         % Max. enzyme transcr. rate         [molecs/min/cell] *
wq = 948.9;         % Max. q-transcr. rate              [molecs/min/cell] *
wp = 0;             % PLACEHOLDER - max. transcr. rate  [molecs/min/ cell]
thetar = 426.9;     % Ribosomal transcr. threshold      [molecs/cell] *
thetax = 4.380;     % Non-ribosomal transcr. threshold  [molecs/cell] *
Kq = 152220;        % q-autoinhibition threshold        [molecs/cell] *
nq = 4;             % q-autoinhibition Hill coeff.      []
M = 1e8;            % Total cell mass                   [aa]

parameters = [s0 ns nr nx gmax Kp vet Ket vem Kem wr we wq wp thetar thetax Kq nq M];

% Define rate constants
dm = 0.1;       % mRNA degradation rate
kb = 1;         % mRNA-ribosome binding rate
ku = 1;         % mRNA-ribosome unbinding rate
rates = [dm kb ku];



%% ODE CONDITIONS

% Initial conditions
r_0 = 10;       % Free ribosomes
et_0 = 0;       % Free transporters
em_0 = 0;       % Free metabolic
q_0 = 0;        % Free housekeeping
p_0 = 0;        % Free PLACEHOLDER
mr_0 = 0;       % mRNA - ribosomal
mt_0 = 0;       % mRNA - transporters
mm_0 = 0;       % mRNA - metabolic
mq_0 = 0;       % mRNA - housekeeping
mp_0 = 0;       % mRNA - PLACEHOLDER
cr_0 = 0;      % Complex - ribosomal
ct_0 = 0;      % Complex - transporters
cm_0 = 0;      % Complex - metabolic
cq_0 = 0;      % Complex - housekeeping
cp_0 = 0;      % Complex - PLACEHOLDER
si_0 = 0;       % Intracellular nutrients
a_0 = 1000;     % Energy

init = [r_0 et_0 em_0 q_0 p_0 mr_0 mt_0 mm_0 mq_0 mp_0 cr_0 ct_0 ...
        cm_0 cq_0 cp_0 si_0 a_0];

% Timespan
t0 = 0;
tf = 1e9;



%% ODE SOLVER

% Call ODE solver
[t,y] = ode15s(@(t,y) Annotate_cellmodel_odes(t, y, rates, parameters), [t0 tf], init);

% Extract variables
r = y(:,1);
et = y(:,2);
em = y(:,3);
q = y(:,4);
p = y(:,5);
mr = y(:,6);
mt = y(:,7);
mm = y(:,8);
mq = y(:,9);
mp = y(:,10);
cr = y(:,11);
ct = y(:,12);
cm = y(:,13);
cq = y(:,14);
cp = y(:,15);
si = y(:,16);
a = y(:,17);

m_all = mr+mt+mm+mq+mp;
c_all = cr+ct+cm+cq+cp;

% Calculate rates from the simulation data
Kgamma = gmax/Kp;
gamma = gmax*a./(Kgamma + a);
ttrate = gamma.*c_all;    
lam = ttrate/M;    
nuimp = et*vet*s0/(Ket + s0);
nucat = em*vem.*si./(Kem + si);