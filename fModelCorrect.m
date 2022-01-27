function [Y,sol,x,t] = fModelCorrect
 global alpha
 k = 54;
 rho = 7800;
 cp = 490;
 alpha = k/(rho*cp);
 m = 0;
 x = linspace(0,0.05,5);
 t = linspace(0,20,200);
 sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
 
 function[c,f,s]=pdex1pde(x,t,u,dudx) % initialising PDE
  c = 1/alpha;
  f = dudx;
  s = 0;
 end

 function u0 = pdex1ic(x)             % Initial Condition
  u0 =20;
 end

 function [pl,ql,pr,qr] =pdex1bc(xl,ul,xr,ur,t)   % Boundary Condition
  pl = ul - 30*sin(t);                 % p = parameter vector or scalar
  ql = 0;
  pr = ur;
  qr = 0;
 end

u = sol(:,:,1) + 0.1*randn(size(sol)); % noisy solution for all x locations
Y = u(:,3,1); % vector of temperature guess values at 3rd x loc 
end

