function ode45_ex_ase166m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: ode45_ex_ase166m.m
%BY: John A. Christian
%CREATED: Feb 23, 2009
%
%This is an example of how to use the ode45 function. It goes along with
%the "Notes on using the ode45 function in MATLAB" PDF I handed out in
%lab on Feb 23, 2009.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear workspace and variables
clear all
close all
clc


%Time will go from 0 to 3 seconds. I want 100 data points.
   tspan = linspace(0,3,100);

%Define initial conditions
   v0 = 20;
   theta0 = 40*pi/180;
   s0 = [0             ;
         0             ;
         v0*cos(theta0);
         v0*sin(theta0) ]

%Call ode45
   options = odeset('RelTol',1e-8, 'AbsTol', 1e-8);
   [t,s] = ode45(@odefun,tspan,s0,options);

%Plot results
   plot(s(:,1),s(:,2)), grid on
   xlabel('x-dimension'), ylabel('y-dimension')
   axis equal

end
   

%% Define the ODEs
function [sdot] = odefun(t,s)

   g0 = 9.81;

   sdot = zeros(4,1);
   
   sdot(1) = s(3);
   sdot(2) = s(4);
   
   sdot(3) = 0;
   sdot(4) = -g0;
   end
   