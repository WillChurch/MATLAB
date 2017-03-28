syms x y t mb lb l mw g Cb Cw Ib Iw T K u

diff(x,t,2) == ((mb*lb + mw*l)*g*sin(x) - K*u - Cb*diff(x,t) + diff(y,t))/(Ib + mw*l^2)
diff(y,t,2) == ((Ib + Iw + mw*l^2)*(T - Cw*diff(y,t))/(Iw*(Ib + mw*l^2))) - (((mb*lb + mw*l)*g*sin(x) - Cb*diff(x,t))/(Ib*mw*l^2))