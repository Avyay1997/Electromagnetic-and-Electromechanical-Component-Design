% Prob F
% To Reproduce the study of ex 6.8A
%
% Written By
% Avyay Sah
% Email: asah@purdue.edu

clc;clear all;close all;

% Given input
ak  = 11.424;                  % material constant (A/m) 
bk  = 0.156;                   % material constant (T^2) 
 c  = 12498;
as  = 0.187;                   % material constant 
bs  = 0.058;                   % material constant (T^2)
cs  = 0.946;                   % material constant 
mu0 = 4e-7*pi;                 % permeability of free space 
f   = 10000;                   % frequency (Hz)  
w   = 2*pi*f;                  % speed (rad/s)
N   = 10000;

T  = 1/f;                      % time period (s)
t  = linspace(0,T,N);          % time (s)
MP = ferrite_catalog(1);       % material characteristics

%Case1
case1.B = 0.4*cos(w*t);                    % Given flux density (T)
case1.dB = -0.4*w*sin(w*t);   
[mu] = muB(MP,case1.B);                    % compute material permeability
case1.Han = case1.B ./ mu;                 % anhysteritic field intensty
case1.Man = case1.B - mu0 * case1.Han;     % anhysterictic magnetization 
case1.k = ak*exp((-case1.Man.^2)/bk);      % pinning site factor
case1.s = cs + as*exp(-(case1.Man.^2)/bs); % energy storage factor
case1.a = - (1./((1+c)*case1.k*mu0)).*...  % model parameter 1
          abs(case1.dB);
case1.b = (((case1.s -(c/(1+c))).*...      % model parameter 2
          case1.Man + ((1-case1.s).*...
          case1.B))./(case1.k*mu0)).*abs(case1.dB);

% initial conditions for Case1
case1.Mirr = ones(1,N)*0.5;              % Irreversible magnetization
case1.h    = zeros(1,N);                 % time step 
case1.ad   = ones(1,N)*0.1;              
case1.bd   = zeros(1,N);

%Case2
case2.B = 0.1*cos(w*t);                    % Given flux density (T)
case2.dB = -0.1*w*sin(w*t);
[mu] = muB(MP,case2.B);                    % compute material permeability
case2.Han = case2.B ./ mu;                 % anhysteritic field intensty
case2.Man = case2.B - mu0 * case2.Han;     % anhysterictic magnetization
case2.k = ak*exp((-case2.Man.^2)/bk);      % pinning site factor
case2.s = cs + as*exp(-(case2.Man.^2)/bs); % energy storage factor
case2.a = - (1./((1+c)*case2.k*mu0)).*...  % model parameter 1
          abs(case2.dB);
case2.b = (((case2.s -(c/(1+c))).*...      % model parameter 2
          case2.Man + ((1-case2.s).*...
          case2.B))./(case2.k*mu0)).*abs(case2.dB);

% initial conditions for Case2
case2.Mirr = ones(1,N)*0.001;              % Irreversible magnetization
case2.h    = zeros(1,N);                   % time step 
case2.ad   = ones(1,N)*0.1;
case2.bd   = zeros(1,N);


% compute irreversible magnetization
for i=1:(N-1)

%---------Case1---------%    
case1.h(i) = t(i+1)-t(i);      

case1.ad(i)= (2 + case1.h(i)*case1.a(i))...
             /(2- case1.h(i)*case1.a(i+1));
         
case1.bd(i)= case1.h(i)*(case1.b(i)+...
             case1.b(i+1))/(2- case1.h(i)...
             *case1.a(i+1));
         
case1.Mirr(i+1) = case1.ad(i)*case1.Mirr(i)...
                 +case1.bd(i);
             
%-------Case2--------%             
case2.h(i) = t(i+1)-t(i);

case2.ad(i)= (2 + case2.h(i)*case2.a(i))...
             /(2- case2.h(i)*case2.a(i+1));
         
case2.bd(i)= case2.h(i)*(case2.b(i)+...
             case2.b(i+1))/(2- case2.h(i)...
             *case2.a(i+1));
         
case2.Mirr(i+1)= case2.ad(i)*case2.Mirr(i)...
                + case2.bd(i);
end

case1.Mrev = (c/(c+1))*(case1.Man...   % compute reversible magnetization
             - case1.Mirr);        
case1.M = case1.Mrev + case1.Mirr;     % compute magnetization
case1.H =(case1.B-case1.M)/mu0;        % compute field intensity  

case2.Mrev = (c/(c+1))*(case2.Man...   % compute reversible magnetization
             - case2.Mirr);
case2.M = case2.Mrev + case2.Mirr;     % copmute magnetization
case2.H =(case2.B-case2.M)/mu0;        % compute field intensity 

alpha = MP.MSE.alpha;
beta  = MP.MSE.beta;
kh    = MP.MSE.kh;
ke    = MP.MSE.ke;

% Compute power loss for Case1
case1.Bmx = max(case1.B);
case1.Bmn = min(case1.B);
case1.DB  = case1.Bmx-case1.Bmn;
case1.feq = 2*trapz(t,case1.dB.^2)...   % equivalent frequency
            /(case1.DB^2*pi^2);
case1.ph  = kh*case1.feq^(alpha-1)...   % hysteresis loss
            *(case1.DB/2)^beta*f;
fprintf('Power loss density for case 1 = %0.00001f W/m^3\n',case1.ph);

%compute power loss for Case 2
case2.Bmx = max(case2.B);
case2.Bmn = min(case2.B);
case2.DB  = case2.Bmx-case2.Bmn;
case2.feq = 2*trapz(t,case2.dB.^2)...   % equivalent frequency
            /(case2.DB^2*pi^2);
case2.ph  = kh*case2.feq^(alpha-1)...   % hysteresis loss
            *(case2.DB/2)^beta*f;
fprintf('Power loss density for case 2 = %0.00001f W/m^3\n',case2.ph);        

% plotting
figure(1)
plot(case1.H,case1.B,case2.H,case2.B)
xlabel('H, A/m');
ylabel('B, T');
title('B vs H');
grid on
legend('Case 1','Case 2','Location','northwest');
