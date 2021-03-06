% interpolacion_vandermonde: calcula los coeficientes del polinomio de
% interpolación a partir de puntos en el plano.
% autor: Sandro León Pastorini.
% año 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% entradas:
%         x vector de coordenadas sobre el eje x para cada punto.
%         y vector de coordenadas sobre el eje y para cada punto.
% así, cada x(i) se corresponde con un y(i)
% salida: p vector de coeficientes representando el polinomio de la interpolación.

function p = interpolacion_vandermonde (x, y)
  n = length(x);
  b = y;
  A = eye(n, n);
  for i = 1:n
    A(i, :) = x .^(i-1);
  endfor
  A = A'
  x = A \ b';
  p = x;
endfunction
