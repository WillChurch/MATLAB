%syms mb mw Ib Iw Cg Cw l lb g Km 

mb = .419;          %kg
mw = .204;
Ib = 3.34*10^-3 ;   %kg*m^2
Iw = .57*10^-3 ;
Cb = 1.02*10^-3 ;   %kg*m^2*s^-2
Cw = 0.05*10^-3 ;
l = .085 ;          %m
lb = .075 ;
g = 9.81 ;          %m*s^-2
Km = 25.1 *10^-3;   %N*m*A^-1

Ts = .02;           %s

%terms
N1 = -Km*Iw*(Ib+mw*l^2);

D3 = Iw*(Ib+mw*l^2)^2;
D2 = (Ib+mw*l^2)*(Cw*(Ib+Iw+mw*l^2)+Iw*Cb);
D1 = Cb*Cw*(2*Iw+Ib+mw*l^2)-Iw*g*(mb*lb+mw*l)*(Ib+mw*l^2);
D0 = -Cw*g*(mb*lb+mw*l)*(Ib+mw*l^2);

s = tf('s');

Gb = minreal((N1*s)/(D3*s^3 + D2*s^2 + D1*s + D0));
Gw = minreal( ((Ib+Iw+mw*l^2)*Km - Iw*(Cb*s + (mb*lb+mw*l)*g)*Gb) ...
        /(Iw*(Ib+mw*l^2)*s + (Ib+Iw+mw*l^2)*Cw));
    
G = [Gb ; Gw];

G.InputName = 'Motor Current';
G.OutputName = ['ThetaB' ; 'OmegaW']

GSS = ss(G)

M = [GSS.B GSS.A*GSS.B (GSS.A^2)*GSS.B (GSS.A^3)*GSS.B (GSS.A^4)*GSS.B (GSS.A^5)*GSS.B];
rank(M);

J = [-2+2*sqrt(3)*i 0 0 0 0 0;
     0 -2-2*sqrt(3)*i 0 0 0 0;
     0 0 -10 0 0 0;
     0 0 0 -10 0 0;
     0 0 0 0 -15 0;
     0 0 0 0 0 -15];
 
 JJ = poly(J);
 
 Phi =  polyvalm(JJ,GSS.A);
 
 K = [ 0 0 0 0 0 1]*inv(M)*Phi