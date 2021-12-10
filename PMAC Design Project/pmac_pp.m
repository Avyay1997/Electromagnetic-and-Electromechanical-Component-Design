% pmac_pp post process a multi-objective optimal design run for a permanent 
%         magnet ac machine. Used in Chapter 9 Section 11 of "Power 
%         Magnetic Devices: A Multi-Objective Design Approach" by 
%         S.D. Sudhoff
%
% Version/Date:
% June 12, 2013
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

% path for libraries-------------------------------------------------------
path(path,'../');
path(path,'Common Routines');

% get the data file--------------------------------------------------------
filename=input('Please enter file name of multiobjective run data: ','s');
load(filename)
clear f size loss

% do the reportplot--------------------------------------------------------
GAP.rp_lvl=1;
distplot(20,fP,1,GAP);
set(20,'Position',[440 378 round(560*0.75) round(420*0.75*0.5)])

% extract information from parameters--------------------------------------
NS=size(bestparameters,2);
for i=1:NS
   parameters=bestparameters(:,i);
   f=pmac_fitness(parameters,D,0);
   size(i)=f.size;
   mass(i)=f.m;
   loss(i)=f.Plagg;
   mss(i)=f.mss;
   mrs(i)=f.mrs;
   mpm(i)=f.mpm;
   mcd(i)=f.mcd;
   Pr1(i)=f.Pr1;
   Pr2(i)=f.Pr2;
   Pr3(i)=f.Pr3;
   Pr4(i)=f.Pr4;
   Pr5(i)=f.Pr5;
   Pr6(i)=f.Pr6;
   Ps1(i)=f.Ps1;
   Ps2(i)=f.Ps2;
   Ps3(i)=f.Ps3;
   Ps4(i)=f.Ps4;
   Ps5(i)=f.Ps5;
   Ps6(i)=f.Ps6;
   Pc1(i)=f.Pc1;
   Pc2(i)=f.Pc2;
   Pc3(i)=f.Pc3;
   Pc4(i)=f.Pc4;
   Pc5(i)=f.Pc5;
   Pc6(i)=f.Pc6;
   Pl1(i)=f.Pl1;
   Pl2(i)=f.Pl2;
   Pl3(i)=f.Pl3;
   Pl4(i)=f.Pl4;
   Pl5(i)=f.Pl5;
   Pl6(i)=f.Pl6;
   J1(i)=f.J1;
   J2(i)=f.J2;
   J3(i)=f.J3;
   J4(i)=f.J4;
   J5(i)=f.J5;
   J6(i)=f.J6;
   P(i)=f.P;
   Is1(i)=f.Is1;
   Is2(i)=f.Is2;
   Is3(i)=f.Is3;
   Is4(i)=f.Is4;
   Is5(i)=f.Is5;
   Is6(i)=f.Is6;
   ac(i)=f.ac;
   st(i)=f.st;
   rt(i)=f.rt;
   ct(i)=f.ct;
   mt(i)=f.mt;
   aslt(i)=f.aslt;
   l(i)=f.l;
   rss(i)=f.rss;
   alphat(i)=f.alphat;
   alphapm(i)=f.alphapm;
   efficiency(i)=f.efficiency;
end

% sort design information by mass------------------------------------------
[size,index]=sort(size,'ascend');
loss=loss(index);
mass=mass(index);
mss=mss(index);
mrs=mrs(index);
mpm=mpm(index);
mcd=mcd(index);
Pr1=Pr1(index);
Pr2=Pr2(index);
Pr3=Pr3(index);
Pr4=Pr4(index);
Pr5=Pr5(index);
Pr6=Pr6(index);
Ps1=Ps1(index);
Ps2=Ps2(index);
Ps3=Ps3(index);
Ps4=Ps4(index);
Ps5=Ps5(index);
Ps6=Ps6(index);
Pc1=Pc1(index);
Pc2=Pc2(index);
Pc3=Pc3(index);
Pc4=Pc4(index);
Pc5=Pc5(index);
Pc6=Pc6(index);
Pl1=Pl1(index);
Pl2=Pl2(index);
Pl3=Pl3(index);
Pl4=Pl4(index);
Pl5=Pl5(index);
Pl6=Pl6(index);
J1=J1(index);
J2=J2(index);
J3=J3(index);
J4=J4(index);
J5=J5(index);
J6=J6(index);
P=P(index);
Is1=Is1(index);
Is2=Is2(index);
Is3=Is3(index);
Is4=Is4(index);
Is5=Is5(index);
Is6=Is6(index);
ac=ac(index);
st=st(index);
rt=rt(index);
ct=ct(index);
mt=mt(index);
aslt=aslt(index);
l=l(index);
rss=rss(index);
alphat=alphat(index);
alphapm=alphapm(index);
efficiency=efficiency(index);
bestparameters=bestparameters(:,index);
index=1:length(index);

% initial set of plots-----------------------------------------------------
figure(4);
plot(size,loss,'gx');
xlabel('Size, sqrt(Kg m^3)');
ylabel('Loss, W');
title('Pareto-Optimal Front');

figure(5)
plot(index,size,'x');
xlabel('Design Number');
ylabel('Size, sqrt(Kg m^3)');
title('Size vs. Design Number');

figure(6)
plot(index,loss,'x');
xlabel('Design Number');
ylabel('Loss, W');
title('Losss vs. Design Number');

figure(7)
plot(mass,J2/mean(J2),'g',mass,Is2/mean(Is2),'m', ...
     mass,ac/mean(ac),'r');
xlabel('Mass, Kg');
legend('Normalized Current Density','Normalized Stator Current', ...
       'Normalized Conductor Area');
title('OP2');
grid on;

figure(8)
plot(mass,Pr2,'g',mass,Ps2,'m',mass,Pc2,'r',mass,Pl2,'b');
xlabel('Mass, kg');
ylabel('Power, W');
legend('Resistive','Semiconductor','Core','Total');
title('OP2');
grid on;

figure(9)
plot(mass,mss,'g',mass,mrs,'m',mass,mcd,'r',mass,mpm,'b');
xlabel('Total Mass, kg');
ylabel('Mass, kg');
legend('Stator Steel','Rotor Steel','Conductor','Magnet');
grid on;

figure(10)
plot(mass,l/mean(l),'g',mass,rss/mean(rss),'m', ...
     mass,alphat,'r',mass,alphapm,'b');
xlabel('Total Mass, kg');
%legend('Normalized Length','Normalized Outer Radius', ...
%       'Tooth Fraction','PM Fraction','Magnet Type');
grid on;

figure(11)
plot(mass,mt,'g',mass,st,'m',mass,rt,'r',mass,ct+0.1,'b');
xlabel('Total Mass, kg');
legend('Magnet','Stator Steel','Rotor Steel','Conductor');
axis([2 10 0 5]);
grid on;


% look at a particular solution--------------------------------------------
dn=input('Enter design number of solution to report on (0 to skip): ');

if (dn>0)

   if (dn>length(index))
      
       disp('No such design');
   
   else

       % show location of design on pareto optimal front
       figure(2);
       plot(size,loss,'bx',size(dn),loss(dn),'ro');
       xlabel('size, sqrt(Kg m^3)');
       ylabel('Loss, W');
       title('Pareto-Optimal Front Loss VS Size');
       legend('All Solutions',['Design ' num2str(dn)]);
       grid on;
      
       % show location of design on pareto optimal front
       figure(200);
       plot(size,efficiency,'bx',size(dn),efficiency(dn),'ro');
       xlabel('size, sqrt(Kg m^3)');
       ylabel('Aggregate Efficiency');
       title('Pareto-Optimal Front Aggregate Efficiency VS Size');
       legend('All Solutions',['Design ' num2str(dn)]);
       grid on;
       
       % give a detailed report on the design
       pmac_fitness(bestparameters(:,dn),D,13);

       GAP.op_list=[1 2];
       reportplot(GAP,GAS,fP);
       
   end 
   
end