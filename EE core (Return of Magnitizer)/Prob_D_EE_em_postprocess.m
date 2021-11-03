%  Prob_D_EE_em_postprocess performs postprocessing of a multi objective 
%                    optimal design of a EI core electromagnet
%                    based on an MEC analysis.  Used in Section 5.4
%                    of "Power Magnetic Devices: A Multi-Objective Design
%                    Approach" by S.D. Sudhoff
%
% Version/Date:
% May 22, 2012
%
% Written by:
% S.D. Sudhoff                               
% Purdue University
% Electrical Engineering Building
% 465 Northwestern Avenue
% West Lafayette, IN 47907-2035
% sudhoff@ecn.purdue.edu
% 765-497-7648
%
% Modified by Avyay Sah
% Email: asah@purdue.edu

% initialize---------------------------------------------------------------
close all
clear all

% get the data file--------------------------------------------------------
filename='Prob_D_EE_em_design_results_2000';
load(filename)


% plot the GA results------------------------------------------------------
GAP.rp_lvl=1;
distplot(20,fP,1,GAP);
reportplot(GAP,GAS,fP);
set(20,'Position',[440 378 560 210])

% find mass & loss & sort by mass -----------------------------------------
NS=size(bi,2);                 % number of solutions
SZ=zeros(1,NS);                % size of electromagnet (m^3)
PL=zeros(1,NS);                % inductor loss
ML=zeros(1,NS);                % Mass of Electromagnet (kg)
VL=zeros(1,NS);                % Volume of Electromagnet
J=zeros(1,NS);                 % current density (A/m^2)
Bg=zeros(1,NS);                % Air-gap Flux density(T)
wc=zeros(1,NS);                % width of center core post (m)
we=zeros(1,NS);                % width of end leg of E core (m)
wb=zeros(1,NS);                % width of base of E core (m)
lc=zeros(1,NS);                % lenght of cores (m)
ds=zeros(1,NS);                % depth of slot (m)
ws=zeros(1,NS);                % width of slot (m)

% extract information from parameters -------------------------------------
for i=1:NS
   parameters=bi(:,i);
   f=Prob_D_EE_em_fit(parameters,D,0);
   SZ(i)=f.SZ;
   VL(i)=f.VL;
   PL(i)=f.PL;
   ML(i)=f.ML;
   Bg(i)=f.Bg;
   J(i)=f.J;
   wc(i)=f.EE.wc;
   we(i)=f.EE.we;
   wb(i)=f.EE.wb;
   lc(i)=f.EE.lc;
   ds(i)=f.EE.ds;
   ws(i)=f.EE.ws;
end

% sort design information bye volume --------------------------------------
[SZ,index]=sort(SZ,'ascend');
PL=PL(index);
VL=VL(index);
ML=ML(index);
Bg=Bg(index);
J=J(index);
wc=wc(index);
we=we(index);
wb=wb(index);
lc=lc(index);
ds=ds(index);
ws=ws(index);
bi=bi(:,index);
index=1:length(index);

% % initial set of plots-----------------------------------------------------
figure(1)
plot(index,SZ,'x');
xlabel('Design Number');
ylabel('Size, sqrt(kg m^3)');
title('Size vs. Design Number');
 
figure(2)
plot(index,PL,'x');
xlabel('Design Number');
ylabel('Loss, W');
title('Losss vs. Design Number');

% look at a particular solution--------------------------------------------
dn=input('Enter design number of solution to report on (0 to skip): ');

if (dn>0)

   if (dn>length(index))
      
       disp('No such design');
   
   else

         
       % show location of design on pareto optimal front
       figure(3);
       plot(SZ,PL,'x-',SZ(dn),PL(dn),'ro');
       xlabel('Size, sqrt(kg m^3)');
       ylabel('Loss, W');
       grid on;
       title('Pareto-Optimal Front');
       legend('All Solutions',['Design ' num2str(dn)]);
       axis([0 10 0 500]);
       
       figure(4)
       plot(SZ,J/1e6);
       title('Size vs Current Density');
       xlabel('Size, sqrt(kg m^3)');
       ylabel('J, A/mm^2');
       grid on;
       
       % give a detailed report on the design
       Prob_D_EE_em_fit(bi(:,dn),D,9);
       
   end 
   
end