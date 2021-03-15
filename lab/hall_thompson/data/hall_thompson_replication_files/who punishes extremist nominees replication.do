/********************************************************************
Name: Stata Replication Code for Who Punishes Extremist Nominees?
Author: Andrew B. Hall and Daniel M. Thompson
Date: December 2017
********************************************************************/

* Set the working directory
cd "~/Dropbox/BaseTurnout/replication file"
set matsize 11000

* Bring in the main analysis dataset
use rd_analysis_hs.dta, clear

* Store the cutoff for the minimum distance between the 
* moderate and extremists necessary to enter the sample
sum absdist, d
gl cutoff = r(p50)


/*
Table 1
*/

* Bring in the main analysis dataset 
use rd_analysis_hs, clear

* column 1
reg vote_general treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* column 2
reg vote_general treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* column 3
reg vote_general treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* column 4
rdrobust vote_general rv if absdist>$cutoff, vce(cluster g)

* column 5
reg victory_general treat rv treat_rv  if abs(rv)<.1 & absdist>$cutoff, cluster(g)

* column 6
reg victory_general treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* column 7
reg victory_general treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* column 8
rdrobust victory_general rv if absdist>$cutoff, vce(cluster g)


/*
Figure 2
*/

* Built from scratch in R


/*
Table 2
*/

* Bring in the main analysis dataset 
use rd_analysis_hs, clear

* column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* column 4
rdrobust turnout_party_share rv if absdist>$cutoff, vce(cluster g)


/*
Table 3
*/

* Set up the Catalist outcome variable
use rd_analysis_hs, clear
gen catop_party_share = catop_share_rep if dem==0
replace catop_party_share = catop_share_dem if dem==1

* column 1
reg catop_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* column 2
reg catop_party_share treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* column 3
reg catop_party_share treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* column 4
rdrobust catop_party_share rv if absdist>$cutoff, vce(cluster g)


/*
Table 4
*/

* Bring in the Hall-Snyder analysis file (row 1)
use rd_analysis_hs, clear

* column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* column 2, row 1
reg turnout_party_share treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* column 3, row 1
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* column 4, row 1
rdrobust turnout_party_share rv if absdist>$cutoff, vce(cluster g)


* Bring in the DIME analysis file (row 2)
use rd_analysis, clear
sum absdist, d
local cutoff_dime = r(p50)

* column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<.1 & absdist>`cutoff_dime', cluster(g)

* column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>`cutoff_dime', cluster(g)

* column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>`cutoff_dime', cluster(g)

* column 4
rdrobust turnout_party_share rv if absdist>`cutoff_dime', vce(cluster g)


* Bring in the Dynamic DIME analysis file (row 3)
use rd_analysis_dyn, clear
sum absdist, d
local cutoff_dyn = r(p50)

* column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<.1 & absdist>`cutoff_dyn', cluster(g)

* column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>`cutoff_dyn', cluster(g)

* column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>`cutoff_dyn', cluster(g)

* column 4
rdrobust turnout_party_share rv if absdist>`cutoff_dyn'


* Bring in the DW-DIME analysis file (row 4)
use rd_analysis_dime, clear
sum absdist, d
local cutoff_dwdime = r(p50)

* column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<.1 & absdist>`cutoff_dwdime', cluster(g)

* column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>`cutoff_dwdime', cluster(g)

* column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>`cutoff_dwdime', cluster(g)

* column 4
rdrobust turnout_party_share rv if absdist>`cutoff_dwdime', vce(cluster g)


/*
Figure 3
*/

* Built from scratch in R


/*
Table 5
*/

* Bring in the main analysis dataset
use rd_analysis_hs, clear

* column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff & turnout_opp_party!=. & !inlist(year,2008,2012), cluster(g)

* column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>$cutoff & !inlist(year,2008,2012), cluster(g)

* column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff & !inlist(year,2008,2012), cluster(g)

* column 4
rdrobust turnout_party_share rv if absdist>$cutoff & turnout_opp_party!=. & !inlist(year,2008,2012), vce(cluster g)

* column 5
reg turnout_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff & turnout_opp_party!=. & inlist(year,2008,2012), cluster(g)

* column 6
reg turnout_party_share treat rv rv2 rv3 if absdist>$cutoff & inlist(year,2008,2012), cluster(g)

* column 7
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff & inlist(year,2008,2012), cluster(g)

* column 8
rdrobust turnout_party_share rv if absdist>$cutoff & turnout_opp_party!=. & inlist(year,2008,2012), vce(cluster g)


/*
Table 6
*/


* Bring in the main analysis dataset
use rd_analysis_hs, clear

* Limit the data to cases where we have no missingness
keep if turnout_party != . & turnout_opp_party != .

* Set up the different components of turnout as outcome variables
gen party_total = rep_count if dem==0
replace party_total = dem_count if dem==1
gen opp_party_total = rep_count if dem==1
replace opp_party_total = dem_count if dem==0
gen voters_party = rep_voters if dem==0
replace voters_party = dem_voters if dem==1
gen voters_opp_party = dem_voters if dem==0
replace voters_opp_party = rep_voters if dem==1 

* row 1, column 1
reg turnout_party treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* row 1, column 2
reg turnout_party treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* row 1, column 3
reg turnout_party treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* row 1, column 4
rdrobust turnout_party rv if absdist>$cutoff

* row 1, column 5
reg turnout_opp_party treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* row 1, column 6
reg turnout_opp_party treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* row 1, column 7
reg turnout_opp_party treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* row 1, column 8
rdrobust turnout_opp_party rv if absdist>$cutoff

* row 2, column 1
reg voters_party treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* row 2, column 2
reg voters_party treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* row 2, column 3
reg voters_party treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* row 2, column 4
rdrobust voters_party rv if absdist>$cutoff

* row 2, column 5
reg voters_opp_party treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* row 2, column 6
reg voters_opp_party treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* row 2, column 7
reg voters_opp_party treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* row 2, column 8
rdrobust voters_opp_party rv if absdist>$cutoff

* row 3, column 1
reg party_total treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* row 3, column 2
reg party_total treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* row 3, column 3
reg party_total treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* row 3, column 4
rdrobust party_total rv if absdist>$cutoff

* row 3, column 5
reg opp_party_total treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* row 3, column 6
reg opp_party_total treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* row 3, column 7
reg opp_party_total treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* row 3, column 8
rdrobust opp_party_total rv if absdist>$cutoff


/*
Figure 4
*/

* Bring in the primary analysis dataset
use rd_analysis_hs, clear

* Create a new outcome representing the difference in the
* same party and opposing party turnout rate
gen diff_in_turnout  = turnout_opp_party - turnout_party

* Loop over many bandwidths and store RD estimate
matrix B = J(10000, 3, .)
local count = 1
quietly forvalues j=0.1(.001).5 {	
	reg diff_in_turnout treat rv rv2 rv3 if absdist > $cutoff & abs(rv) < `j', cluster(g)
	matrix B[`count', 1] = `j'*100
	matrix B[`count', 2] = _b[treat]
	matrix B[`count', 3] = _se[treat]
	local count = `count' + 1
}

* Format the output
svmat B
keep if B1 != .
keep B*
rename B1 bw
rename B2 diff_in_y
rename B3 se

* Construct the confidence bands
gen upper = diff_in_y + 1.96*se
gen lower = diff_in_y - 1.96*se

* Save the output
save "turnout_rate_diff_fig.dta", replace


/*
Table A1
*/

* Bring in the primary analysis dataset
use rd_analysis_hs, clear

* Define a safe district
gen safe = .
replace safe = 1 if lag_pv >= 0.6 | lag_pv <= 0.4
replace safe = 0 if safe != 1
replace safe = . if lag_pv == .

* Print the average value of covariates by full sample and 10% bw 
qui foreach v of varlist dem lag_pv safe open_all pres {
	sum `v' if absdist>$cutoff
	local mean_all = r(mean)
	local n_all = r(N)
	sum `v' if absdist>$cutoff & abs(rv)<0.1
	local mean_sample = r(mean)
	local n_sample = r(N)
	noi di "mean (n) of `v' for all cases: " %4.2f `mean_all' " (" %3.0f `n_all' ")"
	noi di "mean (n) of `v' for local cases: " %4.2f `mean_sample' " (" %3.0f `n_sample' ")"
}


/*
Table A2
*/

* Bring in the primary analysis dataset
use rd_analysis_hs, clear

* Column 1 
reg lag_share treat rv treat_rv  if abs(rv) < .1 & absdist > $cutoff, cluster(g)

* Column 2
reg lag_share treat rv rv2 rv3 if absdist > $cutoff, cluster(g)

* Column 3
reg lag_share treat rv rv2 rv3 rv4 rv5 if absdist > $cutoff, cluster(g)

* Column 4
rdrobust lag_share rv if absdist > $cutoff, vce(cluster g)



/*
Figure A1
*/

* Bring in the main analysis data file
use rd_analysis_hs, clear

* Transform the running variable from share to percent
replace rv = 100 * rv

* Set up a matrix for storing output
matrix B = J(50, 13, .)

* Loop over many possible bandwidths
local count = 1
quietly forvalues j=5(1)50 {

	* Run the simple linear regression and store the estimates
	reg turnout_party_share treat treat_rv rv if absdist > $cutoff & abs(rv) <= `j', cluster(g)
	matrix B[`count', 1] = _b[treat]
	matrix B[`count', 2] = _b[treat] - 1.96*_se[treat]
	matrix B[`count', 3] = _b[treat] + 1.96*_se[treat]
	
	* Run the 3rd-order polynomial regression and store the estimates
	reg turnout_party_share treat rv rv2 rv3 if absdist > $cutoff & abs(rv) <= `j', cluster(g)
	matrix B[`count', 4] = _b[treat]
	matrix B[`count', 5] = _b[treat] - 1.96*_se[treat]
	matrix B[`count', 6] = _b[treat] + 1.96*_se[treat]
	
	* Run the 5th-order polynomial regression and store the estimates
	reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist > $cutoff & abs(rv) <= `j', cluster(g)
	matrix B[`count', 7] = _b[treat]
	matrix B[`count', 8] = _b[treat] - 1.96*_se[treat]
	matrix B[`count', 9] = _b[treat] + 1.96*_se[treat]
	
	* Run the CCT procedure and store the estimates
	rdrobust turnout_party_share rv if absdist > $cutoff, h(`j')
	matrix B[`count', 10] = e(tau_cl)
	matrix B[`count', 11] = e(tau_cl) - 1.96*e(se_tau_cl)
	matrix B[`count', 12] = e(tau_cl) + 1.96*e(se_tau_cl)
	matrix B[`count', 13] = `j'
	
	* Increment the counter
	local count = `count' + 1
}

* Clean up and save the output for the R graphics file
svmat B
keep B*
keep if B1 != .
renvars B1 B2 B3 B4 B5 B6 B7 B8 B9 B10 B11 B12 B13 / ///
	est1 lower1 upper1 est2 lower2 upper2 ///
	est3 lower3 upper3 est4 lower4 upper4 bw
reshape long est lower upper, i(bw) j(modelnum)
gen model = "Local Linear" if modelnum==1
replace model = "3rd Order" if modelnum==2
replace model = "5th Order" if modelnum==3
replace model = "CCT" if modelnum==4
drop modelnum
saveold spec_for_r_hs, replace 


/*
Table A3
*/

* Bring in the primary analysis dataset
use rd_analysis_hs, clear

* Set up the contributions variable to include as a covariate in the regressions
gen total_contrib = cand_dwnom_amount_prim_noninc0 if treat==1
replace total_contrib = cand_dwnom_amount_prim_noninc1 if treat==0
gen log_total_contrib = log(total_contrib)

* Column 1
reg turnout_party_share treat rv treat_rv log_total_contrib ///
	if abs(rv)<0.1 & absdist>$cutoff & turnout_opp_party!=., cluster(g)

* Column 2
reg turnout_party_share treat rv rv2 rv3 log_total_contrib if absdist>$cutoff, cluster(g)

* Column 3
reg turnout_party_share treat rv rv2 rv3 log_total_contrib rv4 rv5 if absdist>$cutoff, cluster(g)

* Column 4
rdrobust turnout_party_share rv if absdist>$cutoff & turnout_opp_party!=., covs(log_total_contrib)



/*
Table A4
*/

* Bring in the state legislature-based scaling analysis data
use rd_analysis_shor_new, clear
sum absdist, d
local shor_cutoff = r(p50)

* Column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>`shor_cutoff', cluster(g)

* Column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>`shor_cutoff', cluster(g)

* Column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>`shor_cutoff', cluster(g)

* Column 4
rdrobust turnout_party_share rv if absdist>`shor_cutoff'


/*
Table A5
*/

* Bring in the analysis dataset with outcome variables that omit leaners
use rd_analysis_no_leaners_hs, clear
sum absdist, d
local no_lean_cutoff = r(p50)

* Column 1
reg turnout_party_share treat rv treat_rv  if abs(rv)<0.1 & absdist>`no_lean_cutoff', cluster(g)

* Column 2
reg turnout_party_share treat rv rv2 rv3 if absdist>`no_lean_cutoff', cluster(g)

* Column 3
reg turnout_party_share treat rv rv2 rv3 rv4 rv5 if absdist>`no_lean_cutoff', cluster(g)

* Column 4
rdrobust turnout_party_share rv if absdist>`no_lean_cutoff', vce(cluster g)


/*
Table A6
*/

* Bring in a dataset of vote choice from the CCES
use rd_analysis_vote_choice_hs, clear
sum absdist, d
local vc_cutoff = r(p50)

* Column 1
reg vc_party treat rv treat_rv  if abs(rv)<0.1 & absdist>`vc_cutoff', cluster(g)

* Column 2
reg vc_party treat rv rv2 rv3 if absdist>`vc_cutoff', cluster(g)

* Column 3
reg vc_party treat rv rv2 rv3 rv4 rv5 if absdist>`vc_cutoff', cluster(g)

* Column 4
rdrobust vc_party rv if absdist>`vc_cutoff', vce(cluster g)


/*
Table A7
*/

* Bring in the main analysis dataset for cases with independents present 
use rd_analysis_hs, clear
keep if turnout_ind!=. & ind_count!= . & ind_voters!=.

* Column 1, Row 1
reg turnout_indep treat rv treat_rv  if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* Column 2, Row 1
reg turnout_indep treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* Column 3, Row 1
reg turnout_indep treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* Column 4, Row 1
rdrobust turnout_indep rv if absdist>$cutoff, vce(cluster g)

* Column 1, Row 2
reg ind_voters treat rv treat_rv if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* Column 2, Row 2
reg ind_voters treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* Column 3, Row 2
reg ind_voters treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* Column 4, Row 2
rdrobust ind_voters rv if absdist>$cutoff, vce(cluster g)

* Column 1, Row 3
reg ind_count treat rv treat_rv if abs(rv)<0.1 & absdist>$cutoff, cluster(g)

* Column 2, Row 3
reg ind_count treat rv rv2 rv3 if absdist>$cutoff, cluster(g)

* Column 3, Row 3
reg ind_count treat rv rv2 rv3 rv4 rv5 if absdist>$cutoff, cluster(g)

* Column 4, Row 3
rdrobust ind_count rv if absdist>$cutoff, vce(cluster g)


/*
Table A8
*/

* Define a matrix and a counter for the matrix row
local row = 1
matrix B = J(10000,9,.)

* Run the placebo check on both sides of the true cutoff separately
qui forval positive=0/1{ 

	* Bring in the data for that side of the cutoff
	use rd_analysis_hs, clear
	if `positive'==0 drop if rv>=0
	if `positive'==1 drop if rv<=0
	
	* Set up a list of placebo cutoff values based on the data and loop over them
	local i = 0
	local min = 40
	local max = _N - 40
	levelsof rv
	qui foreach c in `r(levels)' {
		
		* If the cutoff value leaves too little data on either side, skip it
		local i = `i' + 1
		if !inrange(`i', `min', `max') continue
		
		* Set up the variables for the regressions
		gen rv_c = rv - `c'
		gen treat_c = rv_c > 0
		gen treat_rv_c = treat_c*rv_c
		gen rv_c2 = rv_c^2
		gen rv_c3 = rv_c^3
		gen rv_c4 = rv_c^4
		gen rv_c5 = rv_c^5
		
		* Run the local linear regression
		reg turnout_party_share treat_c rv_c treat_rv_c  if abs(rv_c) < .1 & absdist > $cutoff, cluster(g)
		local b1 = _b[treat_c]
		local se1 = _se[treat_c]
		
		* Run the 3rd order polynomial regression
		reg turnout_party_share treat_c rv_c rv_c2 rv_c3 if absdist > $cutoff, cluster(g)
		local b2 = _b[treat_c]
		local se2 = _se[treat_c]
		
		* Run the 5th order polynomial regression
		reg turnout_party_share treat_c rv_c rv_c2 rv_c3 rv_c4 rv_c5 if absdist > $cutoff, cluster(g)
		local b3 = _b[treat_c]
		local se3 = _se[treat_c]
		
		* Run the rdrobust routine
		rdrobust turnout_party_share rv_c if absdist > $cutoff, vce(cluster g)
		local b4 = e(tau_cl)
		local se4 = e(se_tau_cl)
		
		* Store the output in the matrix
		matrix B[`row', 1] = `c'
		matrix B[`row', 2] = `b1'
		matrix B[`row', 3] = `se1'
		matrix B[`row', 4] = `b2'
		matrix B[`row', 5] = `se2'
		matrix B[`row', 6] = `b2'
		matrix B[`row', 7] = `se3'
		matrix B[`row', 8] = `b4'
		matrix B[`row', 9] = `se4'
		drop rv_c* treat_rv_c treat_c 
		local row = `row' + 1
	}
}

* Format the output
svmat B
keep B*
renvars B1-B9 / c beta1 se1 beta2 se2 beta3 se3 beta4 se4
drop if c==.
reshape long beta se, i(c) j(modelnum)
gen model = "Local Linear" if modelnum==1
replace model = "3rd Order" if modelnum==2
replace model = "5th Order" if modelnum==3
replace model = "CCT" if modelnum==4
gen t = beta/se
gen reject = abs(t)>1.96

* Calculate the average effect and rejection rate for the
* placebo across the four models
collapse beta reject, by(model)


/*
Figure A2
*/

* Bring in the primary analysis dataset
use rd_analysis_hs, clear
keep if absdist>$cutoff

* Run the density plot and McCrary test
DCdensity rv, breakpoint(0) generate(Xj Yj r0 fhat se_fhat) 


/*
Table A9
*/

* Bring in the panel dataset
use panel_analysis, clear
keep if year >= 2006
keep if min > 20

* Standardize the variables
foreach var in avgdonorcfscoreincd lag_pv open inc0 inc1 ///
	numgiverstotal0 numgiverstotal1 dist_btw_cands  {
	sum `var'
	replace `var' = (`var' - r(mean))/r(sd)
}
replace vote_g = 100*vote_g

* Run the regression for the distance from average donor specification
reg vote_g dist_to_donor1 dist_to_donor0, r							 // Column 1
reg vote_g dist_to_donor1 dist_to_donor0 avgd lag_pv open inc0 inc1 /// Column 2
	numgiverstotal0 numgiverstotal1, r

* Run the regression for the canes-wrone specification
reg vote_g caneswrone1 caneswrone0, r						   // Column 3
reg vote_g caneswrone1 caneswrone0 avgd lag_pv open inc0 inc1 /// Column 4
	numgiverstotal0 numgiverstotal1, r
	
* Run the regression for the ansolabehere et al candidate positioning specification
	// note we arbitrarily estimate this for dems, estimate identical for reps
reg vote_g midpoint dist_btw_cands, r			// Column 5
reg vote_g midpoint avgd lag_pv open inc0 inc1 /// Column 6
	numgiverstotal0 numgiverstotal1 dist_btw_cands, r


/*
Table A10
*/

* Bring in the panel data
use data/panel_analysis, clear
keep if min > 20

* Standardize and prepare the variables
foreach var in avgdonorcfscoreincd lag_pv open inc0 inc1 ///
	numgiverstotal0 numgiverstotal1 dist_btw_cands {
	sum `var'
	replace `var' = (`var' - r(mean))/r(sd)
}
replace turnout_share1 = turnout_share1*100
	
* Run the regression for the distance from average donor specification
reg turnout_share1 dist_to_donor1 dist_to_donor0, r				   // Column 1
reg turnout_share1 dist_to_donor1 dist_to_donor0 avgd lag_pv open /// Column 2
	inc0 inc1 numgiverstotal0 numgiverstotal1, r

* Run the regression for the canes-wrone specification
reg turnout_share1 caneswrone1 caneswrone0, r				 // Column 3
reg turnout_share1 caneswrone1 caneswrone0 avgd lag_pv open /// Column 4
	inc0 inc1 numgiverstotal0 numgiverstotal1 , r	

* Run the regression for the ansolabehere et al candidate positioning specification
	// note we arbitrarily estimate this for dems, estimate identical for reps
reg turnout_share1 midpoint dist_btw_cands, r			// Column 5
reg turnout_share1 midpoint avgd lag_pv open inc0 inc1 /// Column 6
	numgiverstotal0 numgiverstotal1 dist_btw_cands, r


/*
Table A11
*/

* Bring in the data on election competitiveness and turnout
use "competitiveness_analysis.dta", clear

* Column 1
areg log_votes margin i.year if midterm==1 & year >= 1980, a(d) cluster(d)

* Column 2
areg log_votes margin i.year if midterm==0 & year >= 1980, a(d) cluster(d)

* Column 3
areg log_votes margin i.year if midterm==1 & year >= 2004, a(d) cluster(d)

* Column 4
areg log_votes margin i.year if midterm==0 & year >= 2004, a(d) cluster(d)

