function pdepe_ex_ase166m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FUNCTION: pdepe_ex_ase166m.m
%BY: John A. Christian
%CREATED: Nov 10, 2008
%
%This is an example of how to use the pdepe function. It goes along with
%the "Notes on using the pdepe function in MATLAB" PDF I handed out in
%lab on November 10, 2008.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Clear workspace and variables
clear all
close all
clc


%Define global variables
global alpha

%Time will go from 0 to 10. I want 100 data points.
   tspan = linspace(0,100,100);

%The x distance goes from 0 to 5 (i.e. 0 <= x <= 5).
%I want 100 data points
   xmesh = linspace(0,1,100);

%Define constants
   m = 0
   alpha = 5e-3

%Call pdepe
   sol = pdepe(m,@pdefun,@icfun,@bcfun,xmesh,tspan);

%Plot results
   surf(xmesh,tspan,sol)
   grid on
   xlabel('x-dimension'), ylabel('time'), zlabel('u')

   
   
end
%% Define the pde function terms
function [c,f,s] = pdefun(x,t,u,DuDx)

   global alpha
   
   c = (1/alpha);
   f = DuDx;
   s = 0;
end

%% Define the initial condition function
function u0 = icfun(x)
   u0 = 1*ones(length(x),1);

end
%% Define the boundary condition function
function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
   
   g = 2-cos(2*pi*t/100);
   
   pl = ul-g;
   ql = 0;
   
   pr = 0;
   qr = 1;
end
