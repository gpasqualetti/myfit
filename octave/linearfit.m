## Copyright (C) 2015 
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {[@var{p}, @var{dp}, @var{chinorm}, @var{r}] =} linearfit (@var{x}, @var{y}, @var{dy}, @var{dx})
## Analytical linear fit.
##
## linearfit tries to fit:
## @tex
## $y = p(1) x + p(2)$
## @end tex
## @ifnottex
## @math{y = p(1) x + p(2)}.
## @end ifnottex
##
## @var{dx} is optional: if supplied the function does a first fit with only @var{dy} and uses the first estimate of @var{p} to propagate @var{dx} on the y. The second fit uses:
## @tex
## $dz = \sqrt{p(1) dx^2 + dy^2}$
## @end tex
## @ifnottex
## @math{dz = sqrt(p(1) dx^2 + dy^2)}.
## @end ifnottex
##
## The return values of @var{p}, @var{dp}, @var{chinorm} and @var{r} are defined as follows.
##
## @table @var
## @item p
## Angular coefficient and intercept obtained by the fit.
##
## @item dp
## Errors on the angular coefficient and the intercept.
##
## @item chinorm
## Normalized chi squared.
## 
## @item r
## Pearson residuals
## @endtable
## @seealso{ols, gls}
## @end deftypefn

## Author:  Giulio Pasqualetti <me@giuliopasqualetti.it>
## Created: 2015-11-16
## Version: 1.0.0

function [p, dp, chinorm, r] = linearfit (x,y,dy,dx)

  if (nargin<3)
    disp("Not enough parameters");
    print_usage();
    return
  end

  if (size(x)!=size(y) || size(x)!=size(dy))
    disp("Size mismatch");
    print_usage();
    return
  end
  
  if nargin==4
    if (size(dx)!=size(x))
      disp("Size mismatch");
      print_usage();
      return
    else
      [p1, dp1, chinorm1, r] = linearfit(x,y,dy);
      dz=sqrt((p1(1)*dx).^2+dy.^2);
      [p, dp, chinorm, r] = linearfit(x,y,dz);
      return
    end
      
  elseif nargin==3
    s0 = sum(1./(dy.^2));
    s1 = sum(x./(dy.^2));
    s2 = sum(x.^2./(dy.^2));
    sy0 = sum(y./(dy.^2));
    sy1 = sum(x.*y./(dy.^2));

    delta = s0*s2 - s1^2;

    p=ones(2,1);
    p(1)=(s0*sy1-sy0*s1)/delta;
    p(2)=(sy0*s2-s1*sy1)/delta;

    dp=ones(2,1);
    dp(1)=sqrt(s0/delta);
    dp(2)=sqrt(s2/delta);

    r = (y-p(1)*x-p(2))./dy;
    chinorm = sum(r.^2)/(length(x)-2);
  end


endfunction
