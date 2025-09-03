

function outvars = demo_1

  [t,y] = ode45(@gugus,[0 20],[2; 0]);

  plot(t,y(:,1),'-o',t,y(:,2),'-o')
  title('Solution of van der Pol Equation (\mu = 1) with ODE45');
  xlabel('Time t');
  ylabel('Solution y');
  legend('y_1','y_2')

##  figure 2
##  plot(y(:,1),y(:,2))

endfunction

function dydt = gugus(t,y)
%VDP1  Evaluate the van der Pol ODEs for mu = 1
%
%   See also ODE113, ODE23, ODE45.

%   Jacek Kierzenka and Lawrence F. Shampine
%   Copyright 1984-2014 The MathWorks, Inc.

  dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];
endfunction
