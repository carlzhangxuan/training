// q-ISLM, Kato(2008)

//Endogenous variables 
var is q  y  r  k  c  g  z  dp mc ;

// Exogenous variables
varexo eG eZ;

// Paramaters
parameters  sigma alpha beta  delta phi lamda  sdG sdZ  omega1 omega2 omega3 theta psi kappa f1 f2 f0 ;
sigma = 1.5;   
alpha = 0.3;   
beta = 0.99;   
delta = 0.025; 
phi = 0.8;
lamda =2;
sdG=0.2;
sdZ=0.1;
omega1 = 0.7;
omega2 = 0.2;
omega3 = 0.1;
theta = 0.97; // =(1-delta)/( alpha*Ystar/kstar + 1-delta )
psi = 0.75 ;  // range 0.01-0.05

kappa = 0.1;
f0 = 0.1;
f1 = 0.5;
f2 = 1.2;

// Model
model(linear);
 
 is = (1/psi)*q +k ;
 q = (1-theta)*(mc(+1) + y(+1) -k(+1)) + theta*q(+1) -r(+1) +g ; 
 y = omega1*c + omega2*is ;
 mc = (1/sigma)*c +(alpha+lamda)*y/(1-alpha) - (1+lamda)*(z+alpha*k)/(1-alpha);
 c =c(+1) - sigma*r(+1) ; 
 r + dp(+1) = f0*(r(-1)+dp) + f1*y + f2*dp(+1) ;  
 k = (1-delta)*k(-1)+delta*is;
 dp = beta*dp(+1) + kappa*mc +(1-beta)*dp(-1) ;

 z = phi*z(-1)+eZ;
 g = phi*g(-1)+eG;
end;

// Initial conditions
initval;
 is = 0;
 y = 0;
 r = 0;
 k = 0;
 c = 0;
 g = 0;
 z = 0;
 
 q  =0;
 mc =0;
 dp =0;
 
 end;

// Shocks
shocks;
 var eG= sdG^2;
 var eZ = sdZ^2;
end;

// Simulation
check;
stoch_simul(irf=40)  y  is  c  k  q  dp z  ;