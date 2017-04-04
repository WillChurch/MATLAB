%syms mb mw Ib Iw Cg Cw l lb g Km 
% 
% rw = .1014/2;                  %wheel radius        %m
% hb = .0063;                    %body height         %m
% wb = .0381;                    %body width          %m
% lb = 9*.0254;                  %body length         %m
% l = 0.0254*5;                  %radius to wheel     %m
% rb = .0950214;                 %radius from pivot to COM   %m
% 
% mb = .162; %.419;              %body mass          %kg
% mm = .110; %.296; %.300;              %motor mass         %kg
% mw = .350;                     %mass of wheel      %kg
% mmw = mw + mm;                 %mass of wheel plus motor %kg  
% 
% Im = 1.35e-5; %5.4*10^-6;         %motor inertia      %kg*m^2  ??????????
% Ib = ((1/12)*mb*(lb^2+wb^2)) + (mb*rb^2); %.00177408 ;       	%kg*m^2             %kg*m^2
% Iw = .5*mw*(rw^2) + Im;  %.00043835 + Im;    %kg*m^2             %kg*m^2
% 
% Cb = 0; %1.02*10^-3 ;       %body friction    

Cw = 0.05*10^-6 ;       %wheel friction

g = 9.81 ;              %m*s^2
%%%%%%%%%%%%WHEEL PARAMETERS
THETA_W = [ .000237 0 0
            0 .000237 0
            0 0 .000469];       %kg*m^2

m_w = .273185;   %kg

r_w = sqrt(THETA_W(1,1)/m_w); %m

%%%%%%%%%%%%FRAME PARAMETERS

THETA_H = [ 0.011058 0.003733 0.004039
            0.003733 0.011096 0.003954
            0.004039 0.003954 0.010378];
       %kg*m^2

m_h = .808223;   %kg

r_h = sqrt(THETA_H/m_w); %m

Km = 25.1e-3; %.45/2 ;            %N*m/A    

M = m_h*r_h + 3* (m_w*r_w);

THETA = THETA_H - (m_h * r_h^2) + 3 * (THETA_W - m_w*r_w^2);

THETA_HAT = THETA - THETA_W;

z3 = zeros(3);

gamma0 = pi/4;         
beta0 = asin(1/sqrt(3));
alpha0 = 0; 

F = [ 0 (sin(gamma0)/cos(beta0)) (cos(gamma0)/cos(beta0))
    0 cos(gamma0) -sin(gamma0)
    1 ((sin(gamma0)*sin(beta0))/cos(beta0)) ((cos(gamma0)*(sin(beta0)))/cos(beta0))
];

DgDphi = [ 0 cos(beta0) 0
        0 sin(beta0)*sin(gamma0) -cos(beta0)*cos(gamma0)
        0 sin(beta0)*cos(gamma0) cos(beta0)*sin(gamma0)
        ];
    
    
    
A3D = [ z3 F z3
        (inv(THETA_HAT)*M*DgDphi) z3 Cw*inv(THETA_HAT)
        -(inv(THETA_HAT)*M*DgDphi) z3 -Cw*(inv(THETA_HAT) + inv(THETA_W))
       ]

B3D = [ z3
        -inv(THETA_HAT)*Km
        (inv(THETA_HAT) + inv(THETA_W))*Km
        ]

C3D = eye(9);

D3D = zeros(9,3);

pen = ss(A3D,B3D,C3D,D3D);

%M = [B3D A3D*B3D (A3D^2)*B3D (A3D^3)*B3D (A3D^4)*B3D (A3D^5)*B3D (A3D^6)*B3D (A3D^7)*B3D (A3D^8)*B3D]
%rank(M)

 J = [(-2+2*sqrt(3)*i) 0 0 0 0 0 0 0 0
      0 (-2-2*sqrt(3)*i) 0 0 0 0 0 0 0
      0 0 0 0 0 0 0 0 0
      0 0 0 0 0 0 0 0 0
      0 0 0 0 -1 0 0 0 0
      0 0 0 0 0 -1 0 0 0
      0 0 0 0 0 0 0 0 0 
      0 0 0 0 0 0 0 0 0
      0 0 0 0 0 0 0 0 0] ;
  
% P3D = [ -10 -15 -8 -5 -5 -1 -2 -3+2*sqrt(3)*i -3-2*sqrt(3)*i]
 
     %-3+2*sqrt(3)*i -3-2*sqrt(3)*i -10
      
     %-5+2*sqrt(3)*i -5-2*sqrt(3)*i -20
      %  -1+2*sqrt(3)*i -1-2*sqrt(3)*i -30]'

%  
% JJ = poly(J);
%  
% Phi =  polyvalm(JJ,A);
%  
% K3D = [0 0 1]*inv(M)*Phi
% 
% AA = A - B*K;

%K3D = place(A3D,B3D,P3D)
K3D = [-5 -7 -5 -1 -1 -1 0 0 0];
open_system('SSModel3D')
