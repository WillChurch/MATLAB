%Let's do some parameter ID based on:
%K*i = J*alpha + B*omega

Km = 25.1e-3; %.45/2 ;            %N*m/A    

Cw55 = Km *(ActualCurrent55./Speed55); 
Cw65 = Km *(ActualCurrent65./Speed65); 
Cw575 = Km *(ActualCurrent575./Speed575); 
Cw75 = Km *(ActualCurrent75./Speed75); 
figure();
subplot(1,2,1);
plot(Time55,Cw55,Time65,Cw65,Time575,Cw575,Time75,Cw75); grid
ylim([-.0000005 .00001])
subplot(1,2,2);
plot(Time55,Speed55,Time65,Speed65,Time575,Speed575,Time75,Speed75);grid


