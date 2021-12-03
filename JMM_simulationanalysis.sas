************************************************
*12/16
*Shared parameter models for simulated data:
*1)Random intercept and slope, slope is shared
*between outcome and exposure.
*2)Random intercept only, intercept is shared
*between outcome and exposure.
************************************************;
options nosource;
ods listing close;
libname ls [insert path];



data simdata;
set ls.jmmsim_coef1_var0;
time=0;
Y=Y0;
A=A0;
theoutcome=Y;
output;

time=1;
Y=Y1;
A=A1;
theoutcome=Y;
output;

time=2;
Y=Y2;
A=.;
theoutcome=Y;
output;

run;

proc sort data=simdata;
by sim id;
run;


ods listing close;

*****************random intercept and slope, only slope is shared;
proc nlmixed data=simdata qpoints=5;
by sim;
parms beta0=1 beta1=-1 beta2=1 beta30=1 
beta31=2 beta5=-1 logsigmasq=0 /*f(Y) parms*/
alpha0=0 alpha1=0.5 alpha2=-0.75 alpha3=1 alpha4=0.35 
alpha5=0 alpha7=1/*f(A) parms*/
logvarb0=0 logvarb1=-1.20 covb0b1=0.27; /*f(b) parms*/

/*normal likelihood for Y*/
mu = beta0 + beta1*X + beta2*(Y0*(time=1) + Y1*(time=2)) +  beta30*A0*(time>0) +  
 beta31*A1*(time>1) + beta5*time + b0 + b1*X*(1 + A0*(time>0) + A1*(time>1));
pY = pdf('NORMAL',Y,mu,exp(logsigmasq));
llhY = log(py);
if llhY=. then llhY=0;


/*Bernoulli likelihood for A*/
logitA = alpha0 + alpha1*X + alpha2*(Y0*(time=0) + Y1*(time=1)) + alpha3*A0*(time=1)  + 
alpha4*X*(Y0*(time=0) + Y1*(time=1)) + alpha5*time + alpha7*b1;
pA = exp(logitA)/(1 + exp(logitA));
llhA = A*log(pA) + (1-A)*log(1-pA);
if llhA=. then llhA=0;

llik = llhY + llhA;

model theoutcome~general(llik);
*lower triangle of vcv matrix in row order;
random b0 b1 ~normal([0,0],[exp(logvarb0),covb0b1,exp(logvarb1)]) subject=id;
*random re ~normal(0,exp(loggre)) subject=id;

estimate 'A0 effect, time 1' beta30;
estimate 'A0 effect, time 2' beta2*beta30 + beta30;
estimate 'A1 effect, time 2' beta31;

ods output ParameterEstimates=ParameterEstimates
           AdditionalEstimates=AdditionalEstimates;
run;





proc sort data=AdditionalEstimates;
by Label;
run;

ods listing;
proc means data=AdditionalEstimates;
by Label;
var estimate standarderror;
run;




*****************random intercept only;
proc nlmixed data=simdata qpoints=5;
by sim;
parms beta0=1 beta1=-1 beta2=1 beta30=1 
beta31=2 beta5=-1 logsigmasq=0 /*f(Y) parms*/
alpha0=0 alpha1=0.5 alpha2=-0.75 alpha3=1 alpha4=0.35 
alpha5=0 alpha6=0/*f(A) parms*/
logvarb0=0 ; /*f(b) parms*/

/*normal likelihood for Y*/
mu = beta0 + beta1*X + beta2*(Y0*(time=1) + Y1*(time=2)) +  beta30*A0*(time>0) +  
 beta31*A1*(time>1) + beta5*time + b0;
pY = pdf('NORMAL',Y,mu,exp(logsigmasq));
llhY = log(py);
if llhY=. then llhY=0;

/*Bernoulli likelihood for A0 and A1*/
logitA = alpha0 + alpha1*X + alpha2*(Y0*(time=0) + Y1*(time=1)) + alpha3*A0*(time=1)  + 
alpha4*X*(Y0*(time=0) + Y1*(time=1)) + alpha5*time + alpha6*b0;
pA = exp(logitA)/(1 + exp(logitA));
llhA = A*log(pA) + (1-A)*log(1-pA);
if llhA=. then llhA=0;

llik = llhY + llhA;


model theoutcome~general(llik);

random b0 ~normal(0,exp(logvarb0)) subject=id;

estimate 'A0 effect, time 1' beta30;
estimate 'A0 effect, time 2' beta2*beta30 + beta30;
estimate 'A1 effect, time 2' beta31;

ods output ParameterEstimates=ParameterEstimates
           AdditionalEstimates=AdditionalEstimates;
run;




proc sort data=AdditionalEstimates;
by Label;
run;

ods listing;
proc means data=AdditionalEstimates;
by Label;
var estimate standarderror;
run;




