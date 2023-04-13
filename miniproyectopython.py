# -*- coding: utf-8 -*-


# seebatt
import math
def seebatt(v):
  c = [0.0, 0.2,  9.0 /  35.0, 16.0 /  63.0, 25.0 /  99.0, 36.0 / 143.0, 49.0 / 195.0, 64.0 / 255.0, 81.0 / 323.0, 100.0 / 399.0, 121.0 / 483.0, 144.0 / 575.0, 169.0 / 675.0, 196.0 / 783.0, 225.0 / 899.0, 256.0 /1023.0, 289.0 /1155.0, 324.0 /1295.0, 361.0 /1443.0, 400.0 /1599.0, 441.0 /1763.0, 484.0 /1935.0]
  sqrtopv = math.sqrt(1.0 + v)
  eta = v / ( 1.0 + sqrtopv )**2
  delold = 1.0
  termold = c[1]   # * eta
  sum1 = termold
  i = 1
  while ((i <= 20) and (abs(termold) > 0.00000001 )):
        del_  = 1.0 / ( 1.0 + c[i+1]*eta*delold )
        term = termold * (del_ - 1.0)
        sum1 = sum1 + term
        i = i + 1
        delold = del_
        termold = term
  return 1.0/ ((1.0/(8.0*(1.0+sqrtopv))) * ( 3.0 + sum1 / ( 1.0+eta*sum1 ) ) );


# test seebatt
if(seebatt(3.76982021737918) - 8.03007718341174 < 10**-12):
  print("test seebatt passed")
else:
  print("test seebatt failed")


# seebattk
def seebattk(v):
  d = [0.0, 1.0/3.0, 4.0/27.0, 8.0/27.0, 2.0 /9.0,22.0 /   81.0,208.0 /  891.0,340.0 / 1287.0,418.0 / 1755.0,598.0/ 2295.0,
                  700.0/ 2907.0,928.0/ 3591.0,1054.0/ 4347.0,1330.0/ 5175.0,1480.0/ 6075.0,1804.0/ 7047.0,1978.0/ 8091.0,2350.0/ 9207.0,
                  2548.0/10395.0,2968.0/11655.0,3190.0/12987.0,658.0/14391.0]
  sum1 = d[1]
  delold = 1.0
  termold = d[1]
  i = 1
  while (True):
        del_  = 1.0 / ( 1.0 + d[i+1]*v*delold )
        term = termold * (del_ - 1.0)
        sum1 = sum1 + term
        i = i + 1
        delold = del_
        termold = term
        if((i<=20) or (abs(termold)>0.000001)):
          break
  return sum1


# test seebattk
if(seebattk(1.31002860149908) - 0.279155336940309 < 10**-12):
  print("test seebattk passed")
else:
  print("test seebattk failed")


# LAMBERTBATTIN
import numpy as np
def LAMBERTBATTIN(ro, r, dm, Dtsec):
  small = 0.000001;
  mu = 3.986004418e14;   # m3/s2
  y1 = 0;
  magr = np.linalg.norm(r);
  magro = np.linalg.norm(ro);
  CosDeltaNu= np.dot(ro,r)/(magro*magr);
  rcrossr = np.cross(ro,r); #linalg.
  magrcrossr = np.linalg.norm(rcrossr);
  if (dm=='pro'):
      SinDeltaNu = magrcrossr/(magro*magr);
  else:
    SinDeltaNu = -magrcrossr/(magro*magr);
  DNu = math.atan2(SinDeltaNu,CosDeltaNu);
  if DNu < 0.0:
    DNu = 2.0*math.pi+DNu;
  RoR   = magr/magro;
  eps   = RoR - 1.0;
  tan2w = 0.25*eps*eps / ( math.sqrt( RoR ) + RoR * ( 2.0 + math.sqrt( RoR ) ) )
  rp    = math.sqrt( magro*magr )*( (math.cos(DNu*0.25))**2 + tan2w );
  if ( DNu < math.pi ):
    L = ( (math.sin(DNu*0.25))**2 + tan2w ) /( (math.sin(DNu*0.25))**2 + tan2w + math.cos( DNu*0.5 ) );
  else:
    L = ( (math.cos(DNu*0.25))**2 + tan2w - math.cos( DNu*0.5 ) ) / ( (math.cos(DNu*0.25))**2 + tan2w ); 
  m    = mu*Dtsec*Dtsec / ( 8.0*rp*rp*rp );
  x    = 10.0;
  xn   = L;
  chord= math.sqrt( magro*magro + magr*magr - 2.0*magro*magr*math.cos( DNu ) );
  s    = ( magro + magr + chord )*0.5;
  lim1 = math.sqrt(m/L);
  Loops= 1;
  while (1):
    x    = xn;    
    tempx= seebatt(x);
    Denom= 1.0 / ( (1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x) ) );
    h1   = ( L+x )**2 * ( 1.0+ 3.0*x + tempx )*Denom;
    h2   = m*( x - L + tempx )*Denom;
    
    # ----------------------- Evaluate CUBIC ------------------
    b = 0.25*27.0*h2 / ((1.0+h1)**3 );
    if (b < -1.0): # reset the initial condition
        xn = 1.0 - 2.0*l;
    else:
      if (y1 > lim1):
        xn = xn * (lim1/y1);
      else:
        u = 0.5*b / ( 1.0 + math.sqrt( 1.0 + b ) )    
        k2 = seebattk(u)
        y = ( ( 1.0+h1 ) / 3.0 )*(2.0 + math.sqrt( 1.0+b )/( 1.0+2.0*u*k2*k2))
        xn= math.sqrt(((1.0-L)*0.5 )**2 + m/(y*y) ) - ( 1.0+L )*0.5;
    Loops = Loops + 1;
    y1=math.sqrt(m/((L+x)*(1.0+x)) );
    if ((abs(xn-x) < small) and (Loops > 30)):
      break

  a=  mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );

  # ------------------ Find Eccentric anomalies -----------------
  # ------------------------ Hyperbolic -------------------------
  if ( a < -small ):
    arg1 = math.sqrt( s / ( -2.0*a ) );
    arg2 = math.sqrt( ( s-chord ) / ( -2.0*a ) );
    # ------- Evaluate f and g functions --------
    AlpH = 2.0 * asinh( arg1 );
    BetH = 2.0 * asinh( arg2 );
    DH   = AlpH - BetH;
    F    = 1.0 - (a/magro)*(1.0 - cosh(DH) );
    GDot = 1.0 - (a/magr) *(1.0 - cosh(DH) );
    G    = Dtsec - math.sqrt(-a*a*a/mu)*(sinh(DH)-DH);
  else:
    # ------------------------ Elliptical ---------------------
    if ( a > small ):
      arg1 = math.sqrt( s / ( 2.0*a ) );
      arg2 = math.sqrt( ( s-chord ) / ( 2.0*a ) );
      Sinv = arg2;
      Cosv = math.sqrt( 1.0 - (magro+magr-chord)/(4.0*a) );
      BetE = 2.0*math.acos(Cosv);
      BetE = 2.0*math.asin(Sinv);
      if ( DNu > math.pi ):
          BetE= -BetE
      Cosv= math.sqrt( 1.0 - s/(2.0*a) );
      Sinv= arg1;
      am  = s*0.5;
      ae  = math.pi;
      be  = 2.0*math.asin( math.sqrt( (s-chord)/s ) );
      tm  = math.sqrt(am*am*am/mu)*(ae - (be-math.sin(be)));
      if ( Dtsec > tm ):
          AlpE= 2.0*math.pi-2.0*math.asin( Sinv )
      else:
          AlpE= 2.0*math.asin( Sinv )
      DE  = AlpE - BetE;
      F   = 1.0 - (a/magro)*(1.0 - math.cos(DE) );
      GDot= 1.0 - (a/magr)* (1.0 - math.cos(DE) );
      G   = Dtsec - math.sqrt(a*a*a/mu)*(DE - math.sin(DE));
    else:
      # --------------------- Parabolic ---------------------
      arg1 = 0.0;
      arg2 = 0.0;
  vo = [0,0,0]
  v = [0,0,0]
  for i in range(3):
    vo[i] = (r[i] - F*ro[i])/G
    v[i] = (GDot*r[i] - ro[i])/G
  return vo,v


# test LAMBERTBATTIN
r1 = [20.0e6, 20.0e6, 0];
r2 = [-20.0e6, 10.0e6, 0];
tof = 1.0 * 86400;
vL1,vL2 = LAMBERTBATTIN(r1, r2, "retro", tof)
for i in range (3):
  print(vL1[i])
for i in range (3):
  print(vL2[i])