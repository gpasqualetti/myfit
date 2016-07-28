function [p, dp, chinorm, r] = linearfit (x,y,dy,dx)

  if (nargin<3)
    disp('Not enough parameters');
    print_usage();
    return
  end

  if (size(x)!=size(y) || size(x)!=size(dy))
    disp('Size mismatch');
    print_usage();
    return
  end
  
  if nargin==4
    if (size(dx)!=size(x))
      disp('Size mismatch');
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
