%syms mb mw Ib Iw Cg Cw l lb g Km 

rw = .1014/2;                  %wheel radius        %m
hb = .0063;                    %body height         %m
wb = .0381;                    %body width          %m
lb = 9*.0254;                  %body length         %m
l = 0.0254*5;                  %radius to wheel     %m
rb = .0950214;                 %radius from pivot to COM   %m

mm = .110; %.296; %.300;              %motor mass         %kg
mw = .273185;                     %mass of wheel      %kg
mbH = (1.62778-mw) + (.5-.475344); %.419;              %body mass          %kg
%The .5 accounts for the mcu etc
mmw = mw + mm;                 %mass of wheel plus motor %kg  

Im = 1.35e-5; %5.4*10^-6;         %motor inertia      %kg*m^2  ??????????
IbH = 0.017274 + (.5-.475344)*rb^2; 
%additional term is for mcu

%((1/12)*mb*(lb^2+wb^2)) + (mb*rb^2); %.00177408 ;       	%kg*m^2             %kg*m^2
Iw = 0.00047 + Im;  %.00043835 + Im;    %kg*m^2             %kg*m^2

Cb = 0; %1.02*10^-3 ;       %body friction      
Cw = 2*10^-6; %0.05*10^-3 ;       %wheel friction

g = 9.81 ;              %m*s^2

Km = 25.1e-3; %.45/2 ;            %N*m/A    

AH = [ 0 1 0 
    ((mbH*rb + mw*l)*g/(IbH+mmw*l^2)) (-Cb/(IbH+mmw*l^2)) (Cw/(IbH+mmw*l^2))
    (-(mbH*rb + mw*l)*g/(IbH+mmw*l^2)) (Cb/(IbH+mmw*l^2)) (-Cw*(IbH+Iw+mmw*l^2)/(Iw*(IbH+mmw*l^2)))];  

BH = [ 0
    (-Km/(IbH+mmw*l^2))
    ((Km*(IbH+Iw+mmw*l^2))/(Iw*(IbH+mmw*l^2)))];

CH = eye(3);

DH = zeros(3,1);

pen = ss(AH,BH,CH,DH);

MH = [BH AH*BH (AH^2)*BH];
rank(MH);

%  J = [-2+2*sqrt(3)*i 0 0;
%       0 -2-2*sqrt(3)*i 0;
%       0 0 -10] 
  
PH = [-2.5+2*sqrt(3)*i -2.5-2*sqrt(3)*i -5]'

%  
% JJ = poly(J);67
%  
% Phi =  polyvalm(JJ,A);
%  
% K = [0 0 1]*inv(M)*Phi
% 
% AA = A - B*K;

KH = place(AH,BH,PH)

QH = [1 0 0
    0 0 0
    0 0 0];
RH = 100000000000000;

[KLQRH, SH, eH] = lqr(AH,BH,QH,RH);
KLQRH
eH
%open_system('SSModelHeavy')
