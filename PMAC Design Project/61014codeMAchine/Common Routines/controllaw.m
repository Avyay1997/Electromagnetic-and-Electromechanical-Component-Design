function [T] = controllaw(Te_des,wr,P,E,vll_max)
%
% Inputs
% Te_des = torque desired (Nm)
% wr     = electrical rotor position (rad/s)
% P      = numeber of poles
% E        = structure of lumped electrical parameters
%  E.Rs    = stator phase resistance (Ohms)        
%  E.Lq    = q-axis inductance (H)
%  E.Ld    = d-axis inductance (H)
%  E.lm    = d-axis flux linkage due to PM (Vs)
% vll_max  = maximum line to line voltage (V)
%
% Internal
% vll_pk   = peak line to line voltage(V)
% a        = temporary variable to solve quadratic equation
% b        = temporary variable to solve quadratic equation
% c        = temporary variable to solve quadratic equation
% Outputs
% T        = structure for current iq and id
%  T.iqsr  = Iq current (A)
%  T.idsr  = Id current (A)
%  T.Dmmint= Discriminant value 
%
% Written by Avyay Sah
% Email: asah@purdue.edu


% iq and wr for the 6 operating points
iqsr = (4*Te_des)/(3*P*E.lm);

% find id for the 6 operating points
idsr    = zeros(1,length(iqsr));
Dmint = zeros(1,length(iqsr));            % Determinant
for i=1:length(iqsr)
    vll_pk = sqrt(3)*sqrt((E.Rs*iqsr(i)+wr(i)*E.Ld*...
             idsr(i)+wr(i)*E.lm)^2+(E.Rs*idsr(i)-wr(i)...
             *E.Lq*iqsr(i))^2);
          if vll_pk > vll_max                       % -ve id to be injected
              Zq_sq = E.Rs^2 + E.Lq^2;
              Zd_sq = E.Rs^2 + E.Ld^2;
              % solve for id using quadratic equation
              a = Zd_sq;
              b = 2*E.Rs*iqsr(i)*wr(i)*(E.Ld-E.Lq)+2*wr(i)^2*E.lm*E.Ld;
              c = iqsr(i)^2*Zq_sq+wr(i)*E.lm*(wr(i)*E.lm+2*E.Rs*iqsr(i))...
                  -(1/3)*vll_max^2;
              Dmint(i) = b^2 - 4*a*c;
              if (Dmint(i)>=0) 
                  idsr(i) = (-b+sqrt(Dmint(i)))/(2*a);
              end
          else
              idsr(i) = 0;
          end
end
T.iqsr = iqsr;
T.idsr = idsr;
T.Dmint = Dmint;
              
              
              
               