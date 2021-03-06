% interpolacion_newton: calcula los coeficientes del polinomio de
% interpolación a partir de puntos en el plano.
% autor: Sandro León Pastorini.
% año 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% entradas:
%         x vector de coordenadas sobre el eje x para cada punto.
%         y vector de coordenadas sobre el eje y para cada punto.
% así, cada x(i) se corresponde con un y(i)
% salida: p vector de coeficientes representando el polinomio de la interpolación.

function p = interpolacion_newton (puntos)
  n = length(puntos);
  b = puntos(:, 2);
  x = puntos(:, 1);
  A = zeros(n, n);
  for i = 2:n
    A();
  endfor
  A = A'
  x = A \ b';
  p = x;
endfunction
