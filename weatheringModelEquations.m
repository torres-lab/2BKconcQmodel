function xprime = f(t,x,A,B,md,mp)
xprime = zeros(2,1);
xprime(1) = (1-x(1)).^md; %Na
xprime(2) = (A.*(1-x(1)).^md) - (A.*B.*heaviside(x(2)-1).*((x(2)-1).^mp));%Si
