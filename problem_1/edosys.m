#!/usr/bin/octave-cli

## Instituto Tecnológico de Costa Rica
## Área Académica de Ingeniería en Computadores
## CE-3102 Análisis Numérico para Ingeniería
## Prof. Pablo Alvarado
## II Semestre 2018
## Examen Final

## PROBLEMA 1

## NOMBRE: Juan Esteban Navarro Camacho
## CARNE: 201236227

1;

global m = 0.1;   ## Masa de la partícula
global b = 0.05;  ## Coeficiente de atenuación
global k = 1;     ## Constante de Hook

##     
## ## Problema 1.1 ##
## ################## 
## Fuerza aplicada en la partícula
global F=@(x,v,t) -k*x -b*v;  ## <<< Ponga aquí su solución

global G=@(x,v,t) v;  ## <<< Ponga aquí su solución
global H=@(x,v,t) F(x,v,t)/m;  ## <<< Ponga aquí su solución



## Resuelva el sistema atenuado masa resorte usando Euler
## tn el último instante de tiempo
## Dt paso temporal
function [t,x]=eulersys(tn,Dt)

  global m b k F;
  
  t=0:Dt:tn; ## Intervalo de simulación

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));
  

  ## Condiciones iniciales
  x(1)=-1;
  ## suponer que velocidad v(1) = 0 inicial igual cero

  

  ## ################## 
  ## ## Problema 1.2 ##´´
  ## ################## 

  ## Resuelva el sistema de ecuaciones con Euler

  for it = 2:size(t)(2)
    v(it) = v(it-1) + (F(x(it-1),v(it-1),t(it-1))/m)*Dt;
    x(it) = x(it-1) + v(it-1)*Dt;
  endfor
endfunction
figure(1,"name","Euler");
hold off;
[t,x]=eulersys(10,0.05);

plot(t,x,"r;\\Delta t=0.05;");

hold on;
[t,x]=eulersys(10,0.01);

plot(t,x,"g;\\Delta t=0.01;");

[t,x]=eulersys(10,0.001);

plot(t,x,"b;\\Delta t=0.001;");

xlabel("t");
ylabel("x(t)");
axis([0,10,-2,2]);
grid on;

## Resuelva el sistema de ecuaciones con Runge-Kutta 4to orden
## tn Último instante de tiempo
## Dt Paso temporal (delta t)
function [t,x] = rksys(tn,Dt)
  global m b k F H G;

  t=0:Dt:tn; ## Intervalo de simulación

  ## Pre-reserve la memoria utilizada.
  x=zeros(size(t));
  v=zeros(size(t));

  ## Condiciones iniciales
  x(1)=-1;

  ## ################## 
  ## ## Problema 1.4 ##
  ## ################## 

  for it = 2:size(t)(2)
    l_1 =  H( x(it-1), v(it-1), t(it-1))*Dt;
    k_1 =  G( x(it-1), v(it-1), t(it-1))*Dt;
    
    l_2 =  H( x(it-1) + k_1/2 , v(it-1) + l_1/2 , t(it)/2 )*Dt;
    k_2 =  G( x(it-1) + k_1/2 , v(it-1) + l_1/2 , t(it)/2 )*Dt;
    
    l_3 =  H( x(it-1) + k_2/2 , v(it-1) + l_2/2 , t(it)/2 )*Dt;
    k_3 =  G( x(it-1) + k_2/2 , v(it-1) + l_2/2 , t(it)/2 )*Dt;
    
    l_4 =  H( x(it-1) + k_3, v(it-1) + l_3, t(it))*Dt;
    k_4 =  G( x(it-1) + k_3, v(it-1) + l_3, t(it))*Dt;
    
    v(it) = v(it-1) + (1/6)* (l_1 + 2*l_2 + 2*l_3 + k_4);
    x(it) = x(it-1) + (1/6)* (k_1 + 2*k_2 + 2*k_3 + k_4);
  endfor

endfunction

figure(2,"name","RK");
hold off;
[t,x]=rksys(10,0.05);
plot(t,x,"r;\\Delta t=0.05;");

hold on;
[t,x]=rksys(10,0.01);
plot(t,x,"m;\\Delta t=0.01;");

[t,x]=rksys(10,0.001);
plot(t,x,"b;\\Delta t=0.001;");

xlabel("t");
ylabel("x(t)");
axis([0,10,-2,2]);
grid on;


