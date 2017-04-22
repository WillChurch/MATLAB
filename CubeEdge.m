%syms mb mw Ib Iw Cg Cw l lb g Km 

rw = .1014/2;                  %wheel radius        %m
hb = .0063;                    %body height         %m
wb = .0381;                    %body width          %m
lb = 9*.0254;                  %body length         %m
l = 0.0254*5;                  %radius to wheel     %m
rb = .0950214;                 %radius from pivot to COM   %m

mm = .110; %.296; %.300;              %motor mass         %kg
mw = .273185;                     %mass of wheel      %kg
mb = .5;%.475344; %.419;              %body mass          %kg
%The .5 accounts for the mcu
mmw = mw + mm;                 %mass of wheel plus motor %kg  

Im = 1.35e-5; %5.4*10^-6;         %motor inertia      %kg*m^2  ??????????
Ib = 0.006611 + (.5-.475344)*rb^2; 
%additional term is for mcu
%((1/12)*mb*(lb^2+wb^2)) + (mb*rb^2); %.00177408 ;       	%kg*m^2             %kg*m^2
Iw = 0.00047 + Im;  %.00043835 + Im;    %kg*m^2             %kg*m^2

Cb = 1.02*10^-3 ;       %body friction      
Cw =  0.05*10^-3 ; 2*10^-6;      %wheel friction

g = 9.81 ;              %m*s^2

Km = 25.1e-3; %.45/2 ;            %N*m/A    

A = [ 0 1 0 
    ((mb*rb + mw*l)*g/(Ib+mmw*l^2)) (-Cb/(Ib+mmw*l^2)) (Cw/(Ib+mmw*l^2))
    (-(mb*rb + mw*l)*g/(Ib+mmw*l^2)) (Cb/(Ib+mmw*l^2)) (-Cw*(Ib+Iw+mmw*l^2)/(Iw*(Ib+mmw*l^2)))];  

B = [ 0
    (-Km/(Ib+mmw*l^2))
    ((Km*(Ib+Iw+mmw*l^2))/(Iw*(Ib+mmw*l^2)))];

C = eye(3);

D = zeros(3,1);

pen = ss(A,B,C,D);

M = [B A*B (A^2)*B];
rank(M);

%  J = [-2+2*sqrt(3)*i 0 0;
%       0 -2-2*sqrt(3)*i 0;
%       0 0 -10] 
  
P = [-2.5+2*sqrt(3)*i -2.5-2*sqrt(3)*i -10]'

%  
% JJ = poly(J);67
%  
% Phi =  polyvalm(JJ,A);
%  
% K = [0 0 1]*inv(M)*Phi
% 
% AA = A - B*K;

K = place(A,B,P)

Q = [1 0 0
    0 1 0
    0 0 1];
R = 2000;

[KLQR, S, e] = lqr(A,B,Q,R);
%KLQR
%e
%open_system('SSModel')
