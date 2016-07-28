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
import numpy

def linearfit(x,y,dy,*dx):
    x = numpy.array(x)
    y = numpy.array(y)
    dy = numpy.array(dy)
    if dx:
        dx = numpy.array(dx)
        p1, dp1, chinorm1, r1 = linearfit(x,y,dy)
        dz = numpy.sqrt((p1[0]*dx)**2 + dy**2)
        p, dp, chinorm, r = linearfit(x,y,dz)
        print("m = %f +- %f\nq = %f +- %f\nchinorm (%d ndof) = %f\nr =" %(p[0],dp[0],p[1],dp[1],len(x)-2,chinorm))
        [print("\t",k) for k in r]
    else:
        s0 = numpy.sum(1/dy**2)
        s1 = numpy.sum(x/dy**2)
        s2 = numpy.sum(x**2/dy**2)
        sy0 = numpy.sum(y/dy**2)
        sy1 = numpy.sum(x*y/dy**2)
        
        delta = s0*s2 - s1**2
        #print(s0,s1,s2,sy0,sy1)

        p = numpy.ones(2)
        p[0] = (s0*sy1 - sy0*s1)/delta
        p[1] = (sy0*s2 - s1*sy1)/delta

        dp = numpy.ones(2)
        dp[0] = numpy.sqrt(s0/delta)
        dp[1] = numpy.sqrt(s2/delta)

        r = (y-p[0]*x-p[1])/dy
        r = r[0]
        chinorm = numpy.sum(r**2)/(len(x)-2)
    return p, dp, chinorm, r
