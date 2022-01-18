function [ t, y ] = backward_euler_fixed ( f, tspan, y0, n )

%*****************************************************************************80
%
%% backward_euler_fixed uses a fixed-point backward Euler method to solve an ODE.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 March 2021
%
%  Author:
%
%    John Burkardt
%
%  Input:
%
%    function handle f: evaluates the right hand side of the ODE.  
%
%    real tspan(2): the starting and ending times.
%
%    real y0(m): the initial conditions. 
%
%    integer n: the number of steps.
%
%  Output:
%
%    real t(n+1,1), y(n+1,m): the solution estimates.
%
  m = length ( y0 );
  t = zeros ( n + 1, 1 );
  y = zeros ( n + 1, m );

  dt = ( tspan(2) - tspan(1) ) / n;

  it_max = 10;

  t(1,1) = tspan(1);
  y(1,:) = y0(:);

  for i = 1 : n
    tp = t(i,1) + dt; 
    yp(1,:) = y(i,:);
    for j = 1 : it_max
      yp(1,:) = y(i,:) + dt * transpose ( f ( tp, yp(1,:) ) );
    end
    t(i+1,1) = tp;
    y(i+1,:) = yp(1,:);
  end

  return
end

