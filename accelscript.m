syms theta theta1 theta2 g at Ts r 

A = (r*(((theta-theta1)/Ts) - ((theta1-theta2)/Ts))/Ts) + g*theta == at ;
pretty(A)

B = simplify( solve(A, theta));
pretty(B)