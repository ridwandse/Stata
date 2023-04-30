



/*                                                             DISCLAIMER
****************************************************************************************************************************************************************
    
	This DO-file is created by Ridwan Ah Sheikh (Delhi School of Economics) for a practice session at  
                      Dept. of Economics - Central University, Kashmir
                           (All the errors and omissions are mine.)
			   
Date: 04/05/2023
Author: Ridwan Ah Sheikh, ridwan@econdse.org
			

			
			
We Randomly generate treatment rollout years uniformly using simulations in STATA.

Generate a complete panel of 300 units observed in 15 periods */

clear all
timer clear
set seed 10
global T = 15
global I = 300   

global pre  5 
global post 8


global ep event_plot
global g0 "default_look"
global g1 xla(-$pre (1) $post) 
global g2 xt("Periods since the event")
global g3 yt("Average causal effect")
global g  $g1 $g2 $g3
global t "together"

set obs `=$I*$T'

* generate unit-id (i) 

gen i = int((_n-1)/$T )+1

* generate calendar-period (t)

gen t = mod((_n-1),$T )+1					
tsset i t

* Randomly generate treatment rollout years (Ei) uniformly 

gen     Ei = 7 if t==1 & i>=1 & i <=35			
replace Ei = 8 if t==1 & i>=36 & i <=70
replace Ei = 9 if t==1 & i>=71 & i <=105
replace Ei = 10 if t==1 & i>=106 & i <=140
replace Ei = 11 if t==1 & i>=141 & i <=175
replace Ei = 12 if t==1 & i>=176 & i <=210
			
bys i (t): replace Ei = Ei[1]

/*generate "relative time" (K), i.e. the number periods since treated 
         (could be missing if never-treated) */
		 
gen K = t-Ei 

* generate a post-treatment "dummy" D=1 for the year unit is first treated and zero otherwise.
								
gen D = K>=0 & Ei!=. 						


* generate treatment effects(tau) and error term epsilon (eps) 

gen tau = 0                                    
replace tau = 20 if  D==1 & i>=1  & i <=35
replace tau = 60 if  D==1 & i>=36 & i <=70
replace tau = 80 if  D==1 & i>=71 & i <=105
replace tau = 100 if D==1 & i>=106 & i <=140
replace tau = 120 if D==1 & i>=141 & i <=175
replace tau = 140 if D==1 & i>=176 & i <=210
gen eps = rnormal()	

* generate outcome variable (Y) with parallel trends and heterogeneous treatment effects	

gen Y = i/100 + tau*D + 3*eps

/* we can alternatively generate the outcome (Y) with unit and time FE's (i) and (t) respectively.

					gen Y = i + 3*t + tau*D + 3*eps 
					
However, FE's play no role since all methods control for them. */

 
 
 
*    For Callaway and Sant'Anna (2021) we ned to generate a variable 'first_treat'
                      

gen first_treat = cond(Ei==., 0, Ei)          // year a unit recieves a treatment



                             

* Gen leads & lags of treatment 
 
cap drop F_*                             //  F_* are pre-treatment leads
cap drop L_*                            //   L_* are post-treatment lags
forval x = 1/14 {  
	gen F_`x' = K == -`x'
}
forval x = 0/5 {
	gen L_`x' = K ==  `x'
}
rename F_1 ref                    // one year before the treatment (reference year)




save CUK_DiD, replace

                 /*           INSTALL THE NECESSARY PACKAGES           
                 You'll need the following commands (packages) to be installed in STATA:
		 
		- did_imputation (Borusyak et al. 2023): available on SSC
		- did_multiplegt (de Chaisemartin and D'Haultfoeuille 2020): available on SSC
	        - csdid and drdid (Callaway and Sant'Anna 2021): available on SSC 
		- bacondecomp  (Andrew Goodman-Bacon 2021): available on SSC
		- ddtiming (Thomas Goldring)                               */
		
ssc install did_imputation, replace
ssc install did_multiplegt, replace
ssc install csdid, replace
ssc install drdid, replace
ssc install reghdfe, replace
ssc install event_plot, replace
ssc install bacondecomp, replace
net install ddtiming, from(https://tgoldring.com/code)
ssc install ftools, replace

********************************************************************************
*                  Dynamic TWFE OLS estimator
*            (standard errors clustered at unit-level)
********************************************************************************
use CUK_DiD, clear
reghdfe Y L_* F_*, a(i t) cluster(i) 

/*  We can use differet level of clustering in the option above, as:
 "cluster(i t)" -> clustering at unit-time level 
 "cluster(t)" -> clustering at time period level */
 
* event-study plot of corresponding TWFE estimator

event_plot, default_look stub_lag(L_#) stub_lead(F_#) together graph_opt(xtitle("Periods since the event") ytitle("OLS coefficients") xlabel(-10(1)08) ///
	title("TWFE OLS"))
	


* (use Rspike in the Graph editor)

********************************************************************************
* Goodman-Bacon decomposition to identify the biases in TWFE model

********************************************************************************

/* we have 6 timing goups, therefore there are (6^2-6 =30) different 2x2 DiD coefficients
                  Late vs Early - (Forbidden) and Early vs Late (Legtimate)  */
		  
use CUK_DiD, clear
bacondecomp Y D, ddetail 

/*  (OOPS ! Note: Panel must be strongly balanced, it does not work with unbalanced data) 


To see the graph more clearly we can do the same decomposition using "ddtiming"
The advantage of [ddtiming] over [bacondecomp] is that it shows the legends more clearly,
we can use Graph editor and denote "Late vs Early" comparisons by hollow circles
                       Graph editor -> Symbol -> Hollow Circle */


use CUK_DiD, clear
ddtiming Y D, i(i) t(t)

********************************************************************************
* Callaway and Sant'Anna (2021) Estimator
********************************************************************************

use CUK_DiD, clear
csdid Y, ivar(i) time(t) gvar(first_treat) notyet
estat simple, estore(cssimp)  
estat group, estore(csgroup)

estat event, estore(csevent) 

event_plot csevent, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") xlabel(-10(1)08) ///
	title("Callaway and Sant'Anna (2021)")) stub_lag(Tp#) stub_lead(Tm#) together
	
/* Tm#: pre-treatment leads
   Tp#: post-treatment lags	
   By dropping "notyet" in the option, it uses never-treated and not-yet treated as as comparison by default
   By writing "never" it only uses never-treated as comparison. 
   (However, if never-treated group is small, it may have inference problems. The standard errors may be large)*/


********************************************************************************
* de Chaisemartin and D'Haultfoeuille (2020) Estimator
********************************************************************************
use CUK_DiD, clear
did_multiplegt Y i t D, robust_dynamic dynamic(08) placebo(10) breps(20) cluster(i)


event_plot e(estimates)#e(variances), default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") ///
	title("de Chaisemartin and D'Haultfoeuille (2020)") xlabel(-10(1)08)) stub_lag(Effect_#) stub_lead(Placebo_#) together
	
/*Alternatively, we can generate event study plot as:

event_plot e(estimates)#e(variances), stub_lag(Effect_#) stub_lead(Placebo_#) $t $g0 graph_opt($g ti("CD 20") name(gCD, replace))

*/


********************************************************************************
* Borusyak et al. (2023) Imputation Estimator
********************************************************************************
use CUK_DiD, clear
did_imputation Y i t Ei, horizons(0/08) pretrend(10) minn(0) 

event_plot, default_look graph_opt(xtitle("Periods since the event") ytitle("Average causal effect") ///
	title("Borusyak et al. (2021) imputation estimator") xlabel(-10(1)08))
	
/*  INTERESTING POINT TO NOTE:

The standard errors are small in imputation estimator, 
because it uses pre-treatment averages as a base-period outcome,rather than period g-1 as in C&S.This has efficiency gain because it averages over more periods.

(Note: BJS & CD don't want too many pre-periods) 

                                            
					                                        
										      								  
								     THANK YOU VERY MUCH FOR FOR INVITING ME.....!

*/


