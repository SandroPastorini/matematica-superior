function disipador()
  %syms t;

  Tamb = 20+273.15; %20°C
  Ti = 46+273.15; %47°C

  n = 10; %nro de particiones;
  a = 108/(164000);
  L = 0.04; %40mm
  T(1) = Ti;
  T(n) = Tamb;
  %a = hc*P/(K*S);
  %K = 164*w/m;
  h = L/(n-1); %aproximación de la derivada (arbitrario o dado en ejercicio)
  hc = 200;

  u = T(x)-T(a);
  u(0) = T(0)-T(a);

  A = coeficientes(n);

  b = independientes(a, h, T);

  X(1:n+1) = solve(A*t-b);

  % T(x) = T(0)*cosh(sqrt(a)*(L-x))/(cosh(sqrt(a)*L);
  %
  % for j=1:n
  %   du(x(j)) = (u(j+1)-u(j))/h
  % end
  %
  % for j=1:n
  %   d2u(x(j)) = (u(j+1)-2*u(j)+u(j-1))/h^2
  % end
endfunction


function M = coeficientes(n)
  nOnes = ones(n+1, 1);
  A = diag((2+a*h^2)*nOnes, 0)+diag(-1*nOnes(1:n), -1)+diag(-1*nOnes(1:n), 1);
  A(1, 1) = 1; A(n+1, n+1) = 1
  %A = -(2+alfa*h^2)*diag(N).+diag(ones(N, 1), 1).+diag(ones(N, 1), -1);
end

function  b = independientes(a, h, T)
  b = a*h^2*T(a)*ones(n+1, 1);
  b(1) = u(0);
  b(n+1) = 0;
endfunction
