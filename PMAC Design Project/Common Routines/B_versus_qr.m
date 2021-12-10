%  FUNCTION:    [dBqr,Bqr,qr]=B_versus_qr(B,phism,thetarm,P)
%  DESCRIPTION: This routine takes an array of flux density values
%               at different teeth/backiron and rotor positions, and uses
%               this to determine the flux density and its derivative
%               at a point versus rotor position.
%  INPUTS:      B        - flux density (# rotor positions
%                          by # stator positions) (T)
%               phism    - mechanical stator positions (rad)
%               thetarm  - mechanical rotor positions (rad)
%               P        - number of poles
%
%  OUTPUTS:     Bqr      - B field at corresponding electrical
%                          rotor positions (T)
%               dBqr     - derivative of B field with respect to electrical
%                          rotor position at corresponding 
%                          rotor positions (T/rad)
%               qr       - electrical rotor positions (rad)
%
%  INTERNAL:    PP       - number of pole pairs
%               qr       - vector of electrical rotor positions (rad)
%               Nphism   - number of mechanical stator positions to
%                          consider
%               Nthetarm - number of mechanical rotor positions
%               phis     - electrical stator positions (rad)
%               i1       - index variable
%               i2       - index variable
%               I        - index for sorting variables
%               A        - matrix for quadratic approximation
%               b        - matrix for quadratic approximation
%               coef     - coefficients for quadratic approximation

function [dBqr,Bqr,qr] = B_versus_qr(B,phism,thetarm,P)

   %Find B field
   PP=P/2;
   phism=phism-phism(1);
   thetar=thetarm'*PP;
   Nphism=length(phism)/PP;
   Nthetarm=length(thetarm);
   N=Nphism*Nthetarm;
   Bqr(N,1)=0;
   qr(N,1)=0;
   for i=1:Nphism
       phis=phism(i)*PP;
       i1=(i-1)*Nthetarm+1;
       i2=i1+Nthetarm-1;
       qr(i1:i2,1)=thetar-phis;
       Bqr(i1:i2,1)=B(:,i)';  
   end
        
   %Map electrical rotor positions to the interval [0,2*pi]
   twopi=2*pi;
   index=find(qr<0);
   qr(index)=qr(index)+twopi;
   index=find(qr>twopi);
   qr(index)=qr(index)-twopi;
   if (max(qr)>2*pi)|(min(qr)<0)
      error('Incorrect electrical position calculation');
   end
    
    
   %Sort the electrical rotor positions and B field vector
   [qr,I]=sort(qr);
   Bqr=Bqr(I);
   
   %Compute derivative of B field 
   %Use quadratic approximation (ax^2+bx+c)
   dBqr=zeros(size(Bqr));
   for i=1:length(qr)
        if i==1
            A=[(qr(end)-2*pi)^2 qr(end)-2*pi 1; qr(i)^2 qr(i) 1; ...
              qr(i+1)^2 qr(i+1) 1];
            b=[Bqr(end); Bqr(i); Bqr(i+1)];
        elseif ((i<length(qr))&&(i>1)) 
            A=[qr(i-1)^2 qr(i-1) 1; qr(i)^2 qr(i) 1; ...
              qr(i+1)^2 qr(i+1) 1];
            b=[Bqr(i-1); Bqr(i); Bqr(i+1)];
        else
            A=[qr(i-1)^2 qr(i-1) 1; qr(i)^2 qr(i) 1;...
              (qr(1)+2*pi)^2 qr(1)+2*pi 1];
            b=[Bqr(i-1); Bqr(i); Bqr(1)];
        end
        coef=inv(A)*b;
        dBqr(i)=2*coef(1)*qr(i)+coef(2);
   end
   
end