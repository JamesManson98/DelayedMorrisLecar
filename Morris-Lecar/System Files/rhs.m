%rhs for finding steady state
function F=rhs(x,par)
F(1) = 1/par(12)*(par(13)-par(3)*(x(1)-par(9))-par(4)*x(2)*(x(1)-par(11))-par(2)*(1/2+1/2*tanh((x(1)-par(5))/par(6)))*(x(1)-par(10))+1/2*par(14)*(1+tanh((x(1)-par(15))/par(16))));
F(2) = par(1)*cosh(1/2*(x(1)-par(7))/par(8))*(1/2+1/2*tanh((x(1)-par(7))/par(8))-x(2));
