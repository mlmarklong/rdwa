*! 1.0.1 Mark C. Long and Jordan Rooklyn 27dec2022
* Fixed typo in output text...replaced "bootrapped" with "bootstrapped".
* Fixed typo in "local graph1" command...replaced "oon" with "on"
* Replaced "LATE" with "RDTE".
* Added option for computing standard errors using Efron's (2014) bias-corrected smoothed standard error method (Equation 7.25).

* version 1.0.2 August 2023
* When referencing results from Pei's rdmse command, responded to Pei's March 2023 changes in that command, which shifted from eclass to rclass.

* version 1.0.3 May 2024
* Fixed one errant line of code to permit graph 3 to run with greater number of samples than lines (replaced [if "`graph3'"!="off" & `samples'<=`graph3_lines' {] with [if "`graph3'"!="off" & `b'<=`graph3_lines' {]

capture program drop rdwa
program define rdwa, eclass
	version 15.1
	syntax varlist(min=2 max=2) [if] [in] [, c(real 0) p_min(integer 0) p_max(integer 4) rbc kernel(string) samples(integer 200) efron2014 details_off graph1 graph2 graph3 graph3_lines(integer 200)]
	local y: word 1 of `varlist'
	local x: word 2 of `varlist'
	if ("`rbc'"=="") local rbc "off"
	else local rbc "on"
	local iskernel = "`kernel'" ==""
	local kernel = cond(`iskernel',"uni", "`kernel'")
	local kernel = substr(strltrim(strlower("`kernel'")),1,3)
	if ("`efron2014'"=="") local efron2014 "off"
	else local efron2014 "on"
	if ("`details_off'"=="") local details "on"
	else local details "off"
	if ("`graph1'"=="") local graph1 "off"
	else local graph1 "on"
	if ("`graph2'"=="") local graph2 "off"
	else local graph2 "on"
	if ("`graph3'"=="") local graph3 "off"
	else local graph3 "on"
	local lines_changed=0
	if `graph3_lines' > `samples' {
		local graph3_lines=`samples'
		local lines_changed=1
	}
	tempfile user_data
	tempvar index bweight weight treated null fig2_x fig2_x_sort fig2_pr_included inv_x_minus_c
	tempname A append_mata y_star y_star_mata t_star t_star_mata y_star_t_star_mata vce cov sum_cov_sq ones_mata y_star_dev_mata t_star_mean_mata t_star_dev_mata z_star_mata cov_mata bc_mata bc
	* STOP THE ANALYSIS IF CUTOFF IS OUTSIDE OF RANGE OF X
	qui sum `x'
	if `c'<r(min) | `c'>r(max){
		display as error "The cutoff must be within the range of the assignment variable (`x')."
		display as error "Change the value of the cutoff before continuing."
	}
	* STOP THE ANALYSIS IF KERNEL TYPE IS MISSPECIFIED
	else if "`kernel'"!="uni" & "`kernel'"!="tri" {
		display as error "Kernel must be uni or tri" 
		display as error "Default is Uniform Weights"
		display as error "- tri - will, instead, use Triangular Weights"
		display as error "Change the entry before continuing"
	}
	* STOP THE ANALYSIS IF p_min IS NEGATIVE
	else if `p_min'<0 {
		display as error "The minimum polynomial order that will be considered must be an integer >=0."
		display as error "Change the value before continuing."
	}
	* STOP THE ANALYSIS IF p_max<p_min
	else if `p_max'<`p_min' {
		display as error "The maximum polynomial order must be >= the polynomial minimum order."
		display as error "Change the value before continuing."
	}
	* STOP THE ANALYSIS IF samples IS NEGATIVE
	else if `samples'<=0 {
		display as error "The samples option must be an integer >0."
		display as error "Change the value before continuing."
	}	
	* STOP THE ANALYSIS IF lines IS NEGATIVE
	else if `graph3_lines'<=0 {
		display as error "The graph3_lines option must be an integer >0."
		display as error "Change the value before continuing."
	}	
	else {
		display " "
		display as text "rdwa began at $S_TIME  $S_DATE"
		display " "
		display as result 	           "{ul:   Notes                                                                        }"
		display as text       			"Per the command's defaults or the user's settings:"
		display as text       			"(a) the treatment is considered to be equal to 1 for observations with `x'>=`c'."
		display as text       			"(b) the algorithm is considering polynomials of order `p_min' to `p_max'."
		if "`rbc'"=="on" {
			display as text			"(c) the algorithm is computing the bias-corrected estimator"
			display as text			"    (with robust standard errors for the CCT & PLCW results)."	
		}
		else {
			display as text       		"(c) the algorithm is computing the conventional estimator"	
			display as text			"    (with conventional standard errors for the CCT & PLCW results)."	
		}
		if ("`kernel'"=="uni") display as text 	"(d) the algorithm is using a uniform (also known as 'rectantular') kernel."	
		else display as text 			"(d) the algorithm is using a triangular (also known as 'edge') kernel."
		display as text       			"(e) the algorithm will generate estimates from `samples' bootstrapped samples."
		if "`graph3'"=="on" {
			display as text                 "(f) a scatterplot and RD lines graph will be produced at the end of the process."
			if `lines_changed'==1 {
				display as text         "(g) the graph3_lines option is set to the user's selected value of the samples option."
			}
		}
		display " "
		qui {
			preserve
			if "`if'"~="" {
				qui keep `if'
			}
			if "`in'"~="" {
				qui keep `in'
			}
			keep `y' `x' 
			drop if `y'==. | `x'==.
			replace `x'=`x'-`c'			// recenter x around cutoff
			foreach v of var `x' `y' {
				local l`v' : variable label `v'
				if `"`l`v''"' == "" {
					label variable `v' `v'
				}			
			}
			local y_label: variable label `y'
			local x_label: variable label `x'
			gen `treated'=`x'>=0
			if `p_max'>0 {
				forvalues p=1/`p_max' {
					tempvar x`p' treated_x`p'
					gen `x`p''=`x'^`p'
					gen `treated_x`p''=`treated'*`x'^`p'
				}		
			}
			sort `x'
			if "`graph2'"!="off" {
				sum `x'
				local x_min=r(min)+`c'
				local x_max=r(max)+`c'
			}
			gen `index'=_n
			save `user_data'
			local N=_N
			count if `x'<0
			local FULL_obs_l=r(N)
			count if `x'>=0
			local FULL_obs_r=r(N)
			local addspaces1 : display %-15.0gc `FULL_obs_l'
			local addspaces1=strlen(strrtrim("`addspaces1'"))
			local addspaces2 : display %-15.0gc `FULL_obs_r'
			local addspaces2=strlen(strrtrim("`addspaces2'"))
			local addspaces=max(`addspaces1',`addspaces2')
			noi display as result "{ul:    Full Sample            }" _continue
			forvalues s=1/`addspaces' {
				if `s'<`addspaces' {
					noi display as result "{ul: }" _continue
				}
				else {
					noi display as result "{ul: }", _newline(0)
				}
			}
			noi display as text   "Number of obs. < cutoff  = " %-15.0gc `FULL_obs_l'
			noi display as text   "Number of obs. >= cutoff = " %-15.0gc `FULL_obs_r'
			noi display " "
			* Obtain Calonico, Cattaneo, and Titiunik's (2014) result assuming a linear specification.
			local CCT_assumed_p=1
			if "`rbc'"=="on" {
				capture rdrobust `y' `x', deriv(0) c(0) p(1) q(2) kernel(`kernel')
				if _rc==0 {
					local CCT_RDTE = e(tau_bc)
					local CCT_RDTE_se = e(se_tau_rb)
					local CCT_best_bw = e(h_l)
					local CCT_obs_l=e(N_h_l)
					local CCT_obs_r=e(N_h_r)
					local df_reg=`CCT_obs_l'+`CCT_obs_r'-2*(1+1)
					local CCT_RDTE_025=`CCT_RDTE'-`CCT_RDTE_se'*invt(`df_reg',0.975)
					local CCT_RDTE_975=`CCT_RDTE'+`CCT_RDTE_se'*invt(`df_reg',0.975)
				}
			}
			else {
				capture rdrobust `y' `x', deriv(0) c(0) p(1) q(2) kernel(`kernel')
				if _rc==0 {
					local CCT_RDTE = e(tau_cl)
					local CCT_RDTE_se = e(se_tau_cl)
					local CCT_best_bw = e(h_l)
					local CCT_obs_l=e(N_h_l)
					local CCT_obs_r=e(N_h_r)
					local df_reg=`CCT_obs_l'+`CCT_obs_r'-2*(1+1)
					local CCT_RDTE_025=`CCT_RDTE'-`CCT_RDTE_se'*invt(`df_reg',0.975)
					local CCT_RDTE_975=`CCT_RDTE'+`CCT_RDTE_se'*invt(`df_reg',0.975)
				}
			}
			* Obtain Pei et al.'s suggested results.
			local lowest_amse=.
			local PLCW_best_p=.
			local PLCW_best_bw=.
			local PLCW_RDTE=.
			local PLCW_RDTE_se=.
			forvalues p=`p_min'/`p_max' {
				local q=`p'+1
				if "`rbc'"=="on" {
					capture rdrobust `y' `x', deriv(0) c(0) p(`p') q(`q') kernel(`kernel')
					if _rc==0 {
						local tau = e(tau_bc)
						local tau_se = e(se_tau_rb)
						local h_l = e(h_l)
						local b_l = e(b_l)
						local obs=e(N_h_l)+e(N_h_r)
						local obs_l=e(N_h_l)
						local obs_r=e(N_h_r)
						capture rdmse `y' `x', c(0) p(`p') deriv(0) kernel(`kernel') h(`h_l') b(`b_l')
						if _rc==0 {
							local amse=r(amse_bc)
							if `amse'<`lowest_amse' {
								local lowest_amse=`amse'
								local PLCW_best_p=`p'
								local PLCW_best_bw=`h_l'
								local PLCW_RDTE=`tau'
								local PLCW_RDTE_se=`tau_se'
								local df_reg=`obs'-2*(1+`PLCW_best_p')
								local PLCW_RDTE_025=`PLCW_RDTE'-`PLCW_RDTE_se'*invt(`df_reg',0.975)
								local PLCW_RDTE_975=`PLCW_RDTE'+`PLCW_RDTE_se'*invt(`df_reg',0.975)
								local PLCW_obs_l=`obs_l'
								local PLCW_obs_r=`obs_r'
							}
						}
					}
				}
				else {
					capture rdrobust `y' `x', deriv(0) c(0) p(`p') q(`q') kernel(`kernel')
					if _rc==0 {
						local tau = e(tau_cl)
						local tau_se = e(se_tau_cl)
						local h_l = e(h_l)
						local b_l = e(b_l)
						local obs=e(N_h_l)+e(N_h_r)
						capture rdmse `y' `x', c(0) p(`p') deriv(0) kernel(`kernel') h(`h_l') b(`b_l')
						if _rc==0 {
							local amse=r(amse_cl)
							if `amse'<`lowest_amse' {
								local lowest_amse=`amse'
								local PLCW_best_p=`p'
								local PLCW_best_bw=`h_l'
								local PLCW_RDTE=`tau'
								local PLCW_RDTE_se=`tau_se'
								local df_reg=`obs'-2*(1+`PLCW_best_p')
								local PLCW_RDTE_025=`PLCW_RDTE'-`PLCW_RDTE_se'*invt(`df_reg',0.975)
								local PLCW_RDTE_975=`PLCW_RDTE'+`PLCW_RDTE_se'*invt(`df_reg',0.975)
							}
						}
					}
				}
			}
			count if (`x'>=0-`PLCW_best_bw' & `x'<0)
			local PLCW_obs_l=r(N)
			count if (`x'>=0 & `x'<=0+`PLCW_best_bw')
			local PLCW_obs_r=r(N)
			if "`graph3'"!="off" {
				if "`kernel'"=="uni" {
					gen `weight'=1 if (`x'>=0-`PLCW_best_bw' & `x'<0)
					replace `weight'=1 if (`x'>=0 & `x'<=0+`PLCW_best_bw')
				}
				else {
					gen `weight'=(1-abs((`x'-0)/`PLCW_best_bw')) if (`x'>=0-`PLCW_best_bw' & `x'<0)
					replace `weight'=(1-abs((`x'-0)/`PLCW_best_bw')) if (`x'>=0 & `x'<=0+`PLCW_best_bw')
				}
				local rhs `treated'
				if `PLCW_best_p'>0 {
					forvalues P=1/`PLCW_best_p' {
						local rhs `rhs' `x`P'' `treated_x`P''
					}
				}
				capture regress `y' `rhs' [weight=`weight']
				local BETA=_b[_cons]
				local PLCW_L_function="y=`BETA'"
				local PLCW_R_function="y=`BETA'"
				local BETA=_b[`treated']
				local PLCW_R_function="`PLCW_R_function'+`BETA'"
				if `PLCW_best_p'>=1 {
					forvalues P=1/`PLCW_best_p' {
						local BETA=_b[`x`P'']
						local PLCW_L_function="`PLCW_L_function'+`BETA'*(x-`c')^`P'"
						local PLCW_R_function="`PLCW_R_function'+`BETA'*(x-`c')^`P'"
						local BETA=_b[`treated_x`P'']
						local PLCW_R_function="`PLCW_R_function'+`BETA'*(x-`c')^`P'"
					}
				}
				sum `x' if `weight'!=.
				local PLCW_L_range_min=r(min)+`c'
				local PLCW_R_range_max=r(max)+`c'
				drop `weight'
			}
			noi display as result "{ul:    Calonico, Cattaneo, & Titiunik's (2014)  results             }"
			noi display as text   "{ul:Poly.Order  Bandwidth  RDTE       se(RDTE)   [95% conf. interval]}"
			noi display as text %-3.0f `CCT_assumed_p' "         " %-9.0g `CCT_best_bw' "  " %-9.0g `CCT_RDTE' "  " %-9.0g `CCT_RDTE_se' "  " %-9.0g `CCT_RDTE_025' "  " %-9.0g `CCT_RDTE_975'
			noi display as text   " "
			noi display as text   "   Number of obs. <  cutoff and w/in bandwidth: " %-15.0gc `CCT_obs_l'
			noi display as text   "   Number of obs. >= cutoff and w/in bandwidth: " %-15.0gc `CCT_obs_r'
			noi display as text   " "
			noi display as result "{ul:    Pei, Lee, Card, & Weber (2021) results                       }"
			noi display as text   "{ul:Poly.Order  Bandwidth  RDTE       se(RDTE)   [95% conf. interval]}"
			noi display as text %-3.0f `PLCW_best_p' "         " %-9.0g `PLCW_best_bw' "  " %-9.0g `PLCW_RDTE' "  " %-9.0g `PLCW_RDTE_se' "  " %-9.0g `PLCW_RDTE_025' "  " %-9.0g `PLCW_RDTE_975'
			noi display as text   " "
			noi display as text   "   Number of obs. <  cutoff and w/in bandwidth: " %-15.0gc `PLCW_obs_l'
			noi display as text   "   Number of obs. >= cutoff and w/in bandwidth: " %-15.0gc `PLCW_obs_r'
			noi display as text   " "
			clear
			* Now, do Long and Rooklyn's method
			matrix `A'=J(`samples',5,.)
			if "`details'"=="on" {
				noi display as result "{ul:   Long and Rooklyn's results for each bootstrapped sample      }"
				noi display as text   "{ul:Sample Poly.Order Bandwidth  RDTE       se(RDTE)    | ave(RDTEs)}"
			}
			else {
				noi display as text "Processing `samples' bootstrapped samples. One dot per sample:"
				noi display as text "{hline 4}{c +}{hline 3} 1 {hline 3}{c +}{hline 3} 2 {hline 3}{c +}{hline 3} 3 {hline 3}{c +}{hline 3} 4 {hline 3}{c +}{hline 3} 5"
			}
			matrix `t_star'=J(`samples',1,.)
			forvalues b=1/`samples' {
				use `user_data'
				if "`efron2014'"=="on" {
					gen `bweight'=.
					bsample, weight(`bweight')
					sort `index'
					putmata `append_mata'=`bweight'
					mata: `y_star_mata'=st_matrix("`y_star'")
					if `b'==1 {
						mata: `y_star_mata'=`append_mata'
					}
					else {
						mata: `y_star_mata'=(`y_star_mata',`append_mata')
					}
					mata: st_matrix("`y_star'", `y_star_mata')
					drop if `bweight'==0
					expand `bweight'
				}
				else {
					bsample
				}
				* Obtain Pei et al.'s suggested results.
				local lowest_amse=.
				local LR_best_p=.
				local LR_best_bw=.
				local LR_RDTE=.
				local LR_RDTE_se=.
				forvalues p=`p_min'/`p_max' {
					local q=`p'+1
					if "`rbc'"=="on" {
						capture rdrobust `y' `x', deriv(0) c(0) p(`p') q(`q') kernel(`kernel')
						if _rc==0 {
							local tau = e(tau_bc)
							local tau_se = e(se_tau_rb)
							local h_l = e(h_l)
							local b_l = e(b_l)
							local obs=e(N_h_l)+e(N_h_r)
							capture rdmse `y' `x', c(0) p(`p') deriv(0) kernel(`kernel') h(`h_l') b(`b_l')
							if _rc==0 {
								local amse=r(amse_bc)
								if `amse'<`lowest_amse' {
									local lowest_amse=`amse'
									local LR_best_p=`p'
									local LR_best_bw=`h_l'
									local LR_RDTE=`tau'
									local LR_RDTE_se=`tau_se'
									local df_reg=`obs'-2*(1+`LR_best_p')
									local LR_LATE_025=`LR_LATE'-`LR_LATE_se'*invt(`df_reg',0.975)
									local LR_LATE_975=`LR_LATE'+`LR_LATE_se'*invt(`df_reg',0.975)
								}
							}
						}
					}
					else {
						capture rdrobust `y' `x', deriv(0) c(0) p(`p') q(`q') kernel(`kernel')
						if _rc==0 {
							local tau = e(tau_cl)
							local tau_se = e(se_tau_cl)
							local h_l = e(h_l)
							local b_l = e(b_l)
							local obs=e(N_h_l)+e(N_h_r)
							capture rdmse `y' `x', c(0) p(`p') deriv(0) kernel(`kernel') h(`h_l') b(`b_l')
							if _rc==0 {
								local amse=r(amse_cl)
								if `amse'<`lowest_amse' {
									local lowest_amse=`amse'
									local LR_best_p=`p'
									local LR_best_bw=`h_l'
									local LR_RDTE=`tau'
									local LR_RDTE_se=`tau_se'
								}
							}
						}
					}
				}
				count if (`x'>=0-`LR_best_bw' & `x'<0)
				local LR_obs_l=r(N)
				count if (`x'>=0 & `x'<=0+`LR_best_bw')
				local LR_obs_r=r(N)
				if "`graph3'"!="off" & `b'<=`graph3_lines' {
					if "`kernel'"=="uni" {
						gen `weight'=1 if (`x'>=0-`LR_best_bw' & `x'<0)
						replace `weight'=1 if (`x'>=0 & `x'<=0+`LR_best_bw')
					}
					else {
						gen `weight'=(1-abs((`x'-0)/`LR_best_bw')) if (`x'>=0-`LR_best_bw' & `x'<0)
						replace `weight'=(1-abs((`x'-0)/`LR_best_bw')) if (`x'>=0 & `x'<=0+`LR_best_bw')
					}
					local rhs `treated'
					if `LR_best_p'>0 {
						forvalues P=1/`LR_best_p' {
							local rhs `rhs' `x`P'' `treated_x`P''
						}
					}
					capture regress `y' `rhs' [weight=`weight']
					local BETA=_b[_cons]
					local LR_L_function_`b'="y=`BETA'"
					local LR_R_function_`b'="y=`BETA'"
					local BETA=_b[`treated']
					local LR_R_function_`b'="`LR_R_function_`b''+`BETA'"
					if `LR_best_p'>=1 {
						forvalues P=1/`LR_best_p' {
							local BETA=_b[`x`P'']
							local LR_L_function_`b'="`LR_L_function_`b''+`BETA'*(x-`c')^`P'"
							local LR_R_function_`b'="`LR_R_function_`b''+`BETA'*(x-`c')^`P'"
							local BETA=_b[`treated_x`P'']
							local LR_R_function_`b'="`LR_R_function_`b''+`BETA'*(x-`c')^`P'"
						}
					}
					sum `x' if `weight'!=.
					local LR_L_range_min_`b'=r(min)+`c'
					local LR_R_range_max_`b'=r(max)+`c'
					drop `weight'
				}
				matrix `t_star'[`b',1]=`LR_RDTE'			
				matrix `A'[`b',1]=`LR_best_p'
				matrix `A'[`b',2]=`LR_best_bw'
				matrix `A'[`b',3]=`LR_RDTE'
				matrix `A'[`b',4]=`LR_obs_l'
				matrix `A'[`b',5]=`LR_obs_r'
				clear
				if "`details'"=="on" {
					svmat `A', names(A)
					sum A3
					local aveRDTE=r(mean)
					clear
					noi display as text   %-6.0f `b' " " %-3.0f `LR_best_p' "        " %-9.0g `LR_best_bw' "  " %-9.0g `LR_RDTE' "  " %-9.0g `LR_RDTE_se' "   | " %-9.0g `aveRDTE'
				}
				else {
					noi display as text "." _continue
					if `b'==int(`b'/50)*50 {
						noi display `b', _newline(0)
					}
				}
			}
			if "`efron2014'"=="on" {
				mata: `t_star_mata'=st_matrix("`t_star'")
				mata: `y_star_mata'=st_matrix("`y_star'")
				mata: `y_star_t_star_mata'=(`y_star_mata'', `t_star_mata')
				mata: st_matrix("`vce'", variance(`y_star_t_star_mata'))
				local Nplus1=`N'+1
				matrix `cov'=`vce'[`Nplus1',1..`N']
				local Bminus1=`samples'-1
				matrix `cov'=`cov'*`Bminus1'/`samples'
				matrix `sum_cov_sq'=`cov'*`cov''
				local sd_tilda_B=(`sum_cov_sq'[1,1])^0.5					// Equation 3.6 in Efron (2014)
				mata: `ones_mata'=J(`N',`samples',1)
				mata: `y_star_dev_mata'=`y_star_mata'-`ones_mata'
				mata: `t_star_mean_mata'=mean(`t_star_mata')
				mata: `ones_mata'=J(`samples',1,1)
				mata: `t_star_mean_mata'=`t_star_mean_mata'*`ones_mata'
				mata: `t_star_dev_mata'=`t_star_mata'-`t_star_mean_mata'
				mata: `ones_mata'=J(`N',1,1)
				mata: `t_star_dev_mata'=`ones_mata'*`t_star_dev_mata''
				mata: `z_star_mata'=`y_star_dev_mata' :* `t_star_dev_mata'
				mata: `cov_mata'=st_matrix("`cov'")
				mata: `ones_mata'=J(1,`samples',1)
				mata: `cov_mata'=`cov_mata''*`ones_mata'
				mata: `bc_mata'=(`z_star_mata'-`cov_mata'):^2
				mata: `bc_mata'=rowsum(colsum(`bc_mata'))/(`samples'^2)
				mata: st_matrix("`bc'",`bc_mata')
				local sd_tilda_B_bc=(`sd_tilda_B'^2-`bc'[1,1])^0.5				// Equation 7.5 in Efron (2014)
			}
			noi display " ", _newline(0)
			noi display as result 		"{ul:   Long and Rooklyn's results                                            }"
			noi display as text 		"{ul:Poly.Order Freq.  ave(BWs)   sd(BWs)    min(BWs)   max(BWs)  | ave(RDTEs)}"
			svmat `A', names(A)
			forvalues p=`p_min'/`p_max' {
				count if A1==`p' & A3!=.
				local LR_poly_`p'=r(N)
				sum A2 if A1==`p' & A3!=.
				local aveBW=r(mean)
				local sdBW=r(sd)
				local minBW=r(min)
				local maxBW=r(max)
				sum A3 if A1==`p'
				local aveRDTE=r(mean)
				if `p'<`p_max' {
					noi display as text   %-3.0f `p' "        " %-6.0f `LR_poly_`p'' " " %-9.0g `aveBW' "  " %-9.0g `sdBW' "  " %-9.0g `minBW' "  " %-9.0g `maxBW' " | " %-9.0g `aveRDTE' "
				}
				else {
					local fmt_p : display %-3.0f `p'
					local fmt_freq : display %-6.0f `LR_poly_`p''
					local fmt_aveBW : display %-9.0g `aveBW'
					local fmt_sdBW : display %-9.0g `sdBW'
					local fmt_minBW : display %-9.0g `minBW'
					local fmt_maxBW : display %-9.0g `maxBW'
					local fmt_aveRDTE : display %-9.0g `aveRDTE'
					local fmt_sdRDTE : display %-9.0g `sdRDTE'
					noi display as text   "{ul:`fmt_p'        `fmt_freq' `fmt_aveBW'  `fmt_sdBW'  `fmt_minBW'  `fmt_maxBW' | `fmt_aveRDTE'}
				}
			}
			count if A3!=.
			local freq=r(N)
			sum A2 if A3!=.
			local aveBW=r(mean)
			local sdBW=r(sd)
			local minBW=r(min)
			local maxBW=r(max)
			sum A3 
			local aveRDTE=r(mean)
			local sdRDTE=r(sd)
			noi display as text   "Overall:   " %-6.0f `freq' " " %-9.0g `aveBW' "  " %-9.0g `sdBW' "  " %-9.0g `minBW' "  " %-9.0g `maxBW' " | " %-9.0g `aveRDTE' "  " %-9.0g `sdLATE'
			noi display " "
			sum A4
			local aveLR_obs_l=r(mean)
			local minLR_obs_l=r(min)
			local maxLR_obs_l=r(max)
			sum A5
			local aveLR_obs_r=r(mean)
			local minLR_obs_r=r(min)
			local maxLR_obs_r=r(max)
			noi display as text   		"      Average number of obs. <  cutoff and w/in bandwidth: " %-15.0gc `aveLR_obs_l'
			noi display as text   		"      Minimum number of obs. <  cutoff and w/in bandwidth: " %-15.0gc `minLR_obs_l'
			noi display as text   		"      Maximum number of obs. <  cutoff and w/in bandwidth: " %-15.0gc `maxLR_obs_l'
			noi display as text   		"      Average number of obs. >= cutoff and w/in bandwidth: " %-15.0gc `aveLR_obs_r'
			noi display as text   		"      Minimum number of obs. >= cutoff and w/in bandwidth: " %-15.0gc `minLR_obs_r'
			noi display as text   		"      Maximum number of obs. >= cutoff and w/in bandwidth: " %-15.0gc `maxLR_obs_r'
			noi display 			" "
			noi display 			"      {ul:From average RDTE across bootrapped samples}"
			noi display as text   		"      Estimated treatment effect: " `aveRDTE'
			noi display 			" "
			noi display 			"      {ul:From distribution of RDTEs across bootrapped samples}"
			_pctile A3, p(2.5 97.5) 
			local p025 = r(r1)
			local p975 = r(r2)
			noi display as text 	        "      Percentile-based 95% confidence interval: " `p025' " to " `p975'
			if `aveRDTE'>0 {
				count if A3<0
				local p1tail=r(N)/`freq'
				noi display as text   	"      Share of RDTEs <0 (Percentile-based 1-tailed p-value): " `p1tail'
			}
			else {
				count if A3>0 & A3!=.
				local p1tail=r(N)/`freq'
				noi display as text   	"      Share of RDTEs >0 (Percentile-based 1-tailed p-value): "  `p1tail'
			}
			noi display as text	        "      Percentile-based 2-tailed p-value: " 2*`p1tail'
			noi display 			" "
			noi display 			"      {ul:From standard deviation of RDTEs across bootrapped samples}"
			noi display as text   		"      Estimated standard efford of the treatment effect: " `sdRDTE'
			noi display as text	        "      Normal-based 2-tailed p-value: " 2*(1-normal(abs(`aveRDTE')/`sdRDTE'))
			noi display as text	        "      Normal-based 95% confidence interval: " `aveRDTE'-invnormal(0.975)*`sdRDTE' " to " `aveRDTE'+invnormal(0.975)*`sdRDTE'
			if "`efron2014'"=="on" {
			noi display 			" "
			noi display 			"      {ul:From Efron's (2014) bias-corrected smoothed standard dev.}"
			noi display as text   		"      Estimated standard efford of the treatment effect: " `sd_tilda_B_bc'
			noi display as text	        "      Normal-based 2-tailed p-value: " 2*(1-normal(abs(`aveRDTE')/`sd_tilda_B_bc'))
			noi display as text	        "      Normal-based 95% confidence interval: " `aveRDTE'-invnormal(0.975)*`sd_tilda_B_bc' " to " `aveRDTE'+invnormal(0.975)*`sd_tilda_B_bc'
			}
			local p025_normal=`aveRDTE'-invnormal(0.975)*`sdRDTE'
			local p975_normal=`aveRDTE'+invnormal(0.975)*`sdRDTE'
			if "`graph1'"!="off" & `samples'>=2 {
				histogram A3, freq bin(40) scheme(s1color) lcolor(black) fcolor(teal) xtitle("Estimated Regression Discontinuity Treatment Effect") xline(`p025', lpattern(dash) lcolor(dknavy)) xline(`p975', lpattern(dash) lcolor(dknavy)) xline(`p025_normal', lpattern(shortdash) lcolor(gray)) xline(`p975_normal', lpattern(shortdash) lcolor(gray)) xline(`aveRDTE', lcolor(maroon)) normal
				local c_date: display %tdCCYY-NN-DD =daily("`c(current_date)'", "DMY")
				local c_time = c(current_time)
				local c_time_date = "`c_date'"+"_" +"`c_time'"
				local time_string = subinstr("`c_time_date'", " ", "_", .)
				local time_string = subinstr("`time_string'", ":", "_", .)
				local time_string = subinstr("`time_string'", "-", "_", .)
				graph save "rdwa_graph1_`time_string'.gph", replace
				graph export "rdwa_graph1_`time_string'.png", replace width(3000)
				graph export "rdwa_graph1_`time_string'.pdf", replace 
			}
			if "`graph2'"!="off" & `samples'>=2 {
				if "`kernel'"=="uni" {
					gen `fig2_x'=A2
					label var `fig2_x' "`x_label'"
					gsort -`fig2_x'
					gen `fig2_pr_included'=_n/`samples'
					expand 2
					replace `fig2_x'=-`fig2_x' if _n>`samples'
					replace `fig2_x'=`fig2_x'+`c'
					sum `fig2_x'
					local bw_min=r(min)
					local bw_max=r(max)
					gsort `fig2_x'
					gen `fig2_x_sort'=_n
					local newobs=`samples'*2+1
					set obs `newobs'
					replace `fig2_x'=`x_min' if `fig2_x'==.
					replace `fig2_x_sort'=-1 if `fig2_x_sort'==.
					local newobs=`newobs'+1
					set obs `newobs'
					replace `fig2_x'=`x_max' if `fig2_x'==.
					replace `fig2_x_sort'=`newobs' if `fig2_x_sort'==.
					local newobs=`newobs'+1
					set obs `newobs'
					replace `fig2_x'=`bw_min' if `fig2_x'==.
					replace `fig2_x_sort'=0 if `fig2_x_sort'==.
					local newobs=`newobs'+1
					set obs `newobs'
					replace `fig2_x'=`bw_max' if `fig2_x'==.
					replace `fig2_x_sort'=2*`samples'+1 if `fig2_x_sort'==.
					replace `fig2_pr_included'=0 if `fig2_pr_included'==.
					sort `fig2_x_sort'
				}
				else {
					gen `fig2_x'=A2
					label var `fig2_x' "`x_label'"
					keep `fig2_x'
					sort `fig2_x'
					replace `fig2_x'=`fig2_x'+`c'
					gen `inv_x_minus_c'=1/(`fig2_x'-`c')
					sum `inv_x_minus_c'
					gen `fig2_pr_included'=1-r(mean)*(`fig2_x'-`c') if _n==1
					qui forvalues b=2/`samples' {
						sum `inv_x_minus_c' if _n>=`b'
						replace `fig2_pr_included'=r(N)/`samples' - (r(N)*r(mean))*(`fig2_x'-`c')/`samples' if _n==`b'
					}
					expand 2
					replace `fig2_x'=`c'-(`fig2_x'-`c') if _n>`samples'
					local newobs=`samples'*2+1
					set obs `newobs'
					replace `fig2_x'=`c' if `fig2_x'==.
					replace `fig2_pr_included'=1 if `fig2_pr_included'==.
					local newobs=`newobs'+1
					set obs `newobs'
					replace `fig2_x'=`x_min' if `fig2_x'==.
					local newobs=`newobs'+1
					set obs `newobs'
					replace `fig2_x'=`x_max' if `fig2_x'==.
					replace `fig2_pr_included'=0 if `fig2_pr_included'==.
					sort `fig2_x'
				}
				merge 1:1 _n using `user_data'
				gen `null'=.
				twoway (scatter `null' `x') (line `fig2_pr_included' `fig2_x', lcolor(dknavy) legend(off) scheme(s1color) ytitle("Effective Weight Placed on Observations") ylabel(0 "0%" 0.5 "50%" 1 "100%", angle(0)))
				local c_date: display %tdCCYY-NN-DD =daily("`c(current_date)'", "DMY")
				local c_time = c(current_time)
				local c_time_date = "`c_date'"+"_" +"`c_time'"
				local time_string = subinstr("`c_time_date'", " ", "_", .)
				local time_string = subinstr("`time_string'", ":", "_", .)
				local time_string = subinstr("`time_string'", "-", "_", .)
				graph save "rdwa_graph2_`time_string'.gph", replace
				graph export "rdwa_graph2_`time_string'.png", replace width(3000)
				graph export "rdwa_graph2_`time_string'.pdf", replace 
			}
			clear

			ereturn clear
			ereturn scalar FULL_obs_l = `FULL_obs_l'
			ereturn scalar FULL_obs_r = `FULL_obs_r'
			
			ereturn scalar CCT_P = `CCT_assumed_p'
			ereturn scalar CCT_BW = `CCT_best_bw'
			ereturn scalar CCT_RDTE = `CCT_RDTE'
			ereturn scalar CCT_RDTE_se = `CCT_RDTE_se'
			ereturn scalar CCT_95ci_lower = `CCT_RDTE_025'
			ereturn scalar CCT_95ci_upper = `CCT_RDTE_975'
			ereturn scalar CCT_obs_l = `CCT_obs_l'
			ereturn scalar CCT_obs_r = `CCT_obs_r'
			
			ereturn scalar PLCW_P = `PLCW_best_p'
			ereturn scalar PLCW_BW = `PLCW_best_bw'
			ereturn scalar PLCW_RDTE = `PLCW_RDTE'
			ereturn scalar PLCW_RDTE_se = `PLCW_RDTE_se'
			ereturn scalar PLCW_95ci_lower = `PLCW_RDTE_025'
			ereturn scalar PLCW_95ci_upper = `PLCW_RDTE_975'
			ereturn scalar PLCW_obs_l = `PLCW_obs_l'
			ereturn scalar PLCW_obs_r = `PLCW_obs_r'
			
			forvalues p=`p_min'/`p_max' {
				ereturn scalar LR_P`p'_freq = `LR_poly_`p''
			}
			ereturn scalar LR_ave_BWs = `aveBW'
			ereturn scalar LR_sd_BWs = `sdBW'
			ereturn scalar LR_min_BWs = `minBW'
			ereturn scalar LR_max_BWs = `maxBW'

			ereturn scalar LR_ave_obs_l = `aveLR_obs_l'
			ereturn scalar LR_ave_obs_r = `aveLR_obs_r'
			ereturn scalar LR_min_obs_l = `minLR_obs_l'
			ereturn scalar LR_min_obs_r = `minLR_obs_r'
			ereturn scalar LR_max_obs_l = `maxLR_obs_l'
			ereturn scalar LR_max_obs_r = `maxLR_obs_r'

			ereturn scalar LR_ave_RDTEs = `aveRDTE'

			ereturn scalar LR_pct_95ci_lower = `p025'
			ereturn scalar LR_pct_95ci_upper = `p975'
			ereturn scalar LR_pct_p_1tail = `p1tail'
			ereturn scalar LR_pct_p_2tail = 2*`p1tail'

			ereturn scalar LR_sd_RDTEs = `sdRDTE'
			ereturn scalar LR_normal_p_2tail = 2*(1-normal(abs(`aveRDTE')/`sdRDTE'))
			ereturn scalar LR_normal_95ci_lower = `p025_normal'
			ereturn scalar LR_normal_95ci_upper = `p975_normal'
			
			if "`efron2014'"=="on" {
				ereturn scalar LR_efron2014_sd = `sd_tilda_B_bc'
				ereturn scalar LR_efron2014_normal_p_2tail = 2*(1-normal(abs(`aveRDTE')/`sd_tilda_B_bc'))
				ereturn scalar LR_efron2014_normal_95ci_lower = `aveRDTE'-invnormal(0.975)*`sd_tilda_B_bc'
				ereturn scalar LR_efron2014_normal_95ci_upper = `aveRDTE'+invnormal(0.975)*`sd_tilda_B_bc'
			}
			if "`graph3'"!="off" {
				noi display " "
				noi display as result "{ul:    Working on creating the 3rd graph                }"
				use `user_data'
				replace `x'=`x'+`c'			// reset x to original units
				local mcolorpct=min(30,max(1,500/(`N'^0.5)))
				local lcolorpct=min(50,max(1,int(1000/`graph3_lines')))
				local mult=min(1,max(0.07,100/(`N'^0.5)))
				local msize "*`mult'"
				gen byte `null'=.
				if `graph3_lines'==`samples' {
					local text1="Long and Rooklyn's"
					local text2="Specifications"
				}
				else {
					local text1="`graph3_lines' Long and Rooklyn"
					local text2="Specifications (of `samples')"
				}
				twoway 	(scatter `null' `x', msymbol(none)) (scatter `null' `x' if `x'<0, mlwidth(none) mlcolor(gray%30) mfcolor(gray%30) msize(small)) (line `null' `x', lcolor(purple)) (line `null' `x', lcolor(red%50)) (scatter `null' `x', msymbol(none)) (scatter `null' `x' if `x'>=0, mlwidth(none) mlcolor(gray%30) mfcolor(gray%30) msize(small)) (line `null' `x', lcolor(midgreen)) (line `null' `x', lcolor(blue%50) lwidth(medium) scheme(s1color) xsize(6.5) ysize(7.5) ytitle("`y_label'" " ") xtitle("`x_label'" " ") legend(textfirst bmargin(zero) size(vsmall) col(2) colfirst order(1 "{bf:LEFT OF THRESHOLD}" 2 "Observed Data" 3 "Pei, Lee, Card, and" "Weber's Specification" 4 "`text1'" "`text2'" 5 "{bf:RIGHT OF THRESHOLD}" 6 "Observed Data" 7 "Pei, Lee, Card, and" "Weber's Specification" 8 "`text1'" "`text2'")))
				addplot: scatter `y' `x', mlwidth(none) mfcolor(gray%`mcolorpct') msize(`msize') legend(textfirst bmargin(zero) size(vsmall) col(2) colfirst order(1 "{bf:LEFT OF THRESHOLD}" 2 "Observed Data" 3 "Pei, Lee, Card, and" "Weber's Specification" 4 "`text1'" "`text2'" 5 "{bf:RIGHT OF THRESHOLD}" 6 "Observed Data" 7 "Pei, Lee, Card, and" "Weber's Specification" 8 "`text1'" "`text2'"))
				noi display "Adding `graph3_lines' regression lines from bootstrapped samples"
				noi display as text "{hline 4}{c +}{hline 3} 1 {hline 3}{c +}{hline 3} 2 {hline 3}{c +}{hline 3} 3 {hline 3}{c +}{hline 3} 4 {hline 3}{c +}{hline 3} 5"
				local missingline=0
				forvalues b=1/`graph3_lines' {
					if `A'[`b',1]!=. {
						addplot, nodraw: function `LR_L_function_`b'', range(`LR_L_range_min_`b'' `c') lcolor(red%`lcolorpct') legend(textfirst bmargin(zero) size(vsmall) col(2) colfirst order(1 "{bf:LEFT OF THRESHOLD}" 2 "Observed Data" 3 "Pei, Lee, Card, and" "Weber's Specification" 4 "`text1'" "`text2'" 5 "{bf:RIGHT OF THRESHOLD}" 6 "Observed Data" 7 "Pei, Lee, Card, and" "Weber's Specification" 8 "`text1'" "`text2'"))
						addplot, nodraw: function `LR_R_function_`b'', range(`c' `LR_R_range_max_`b'') lcolor(blue%`lcolorpct') legend(textfirst bmargin(zero) size(vsmall) col(2) colfirst order(1 "{bf:LEFT OF THRESHOLD}" 2 "Observed Data" 3 "Pei, Lee, Card, and" "Weber's Specification" 4 "`text1'" "`text2'" 5 "{bf:RIGHT OF THRESHOLD}" 6 "Observed Data" 7 "Pei, Lee, Card, and" "Weber's Specification" 8 "`text1'" "`text2'"))
						noi display as text "." _continue
					}
					else {
						noi display as text "*" _continue
						local missingline=1
					}
					if `b'-int(`b'/50)*50 == 0 {
						noi display `b', _newline(0)
					}
				}
				if `missingline'==1 {
					noi display as text "* denotes a missing line as the Pei et al. algorithm failed to identify any specifications for this bootstrapped sample"
				}
				addplot, nodraw: function `PLCW_L_function', range(`PLCW_L_range_min' `c') lcolor(purple) legend(textfirst bmargin(zero) size(vsmall) col(2) colfirst order(1 "{bf:LEFT OF THRESHOLD}" 2 "Observed Data" 3 "Pei, Lee, Card, and" "Weber's Specification" 4 "`text1'" "`text2'" 5 "{bf:RIGHT OF THRESHOLD}" 6 "Observed Data" 7 "Pei, Lee, Card, and" "Weber's Specification" 8 "`text1'" "`text2'"))
				addplot: function `PLCW_R_function', range(`c' `PLCW_R_range_max') lcolor(midgreen) legend(textfirst bmargin(zero) size(vsmall) col(2) colfirst order(1 "{bf:LEFT OF THRESHOLD}" 2 "Observed Data" 3 "Pei, Lee, Card, and" "Weber's Specification" 4 "`text1'" "`text2'" 5 "{bf:RIGHT OF THRESHOLD}" 6 "Observed Data" 7 "Pei, Lee, Card, and" "Weber's Specification" 8 "`text1'" "`text2'"))
				noi display, _newline(0)
				local c_date: display %tdCCYY-NN-DD =daily("`c(current_date)'", "DMY")
				local c_time = c(current_time)
				local c_time_date = "`c_date'"+"_" +"`c_time'"
				local time_string = subinstr("`c_time_date'", " ", "_", .)
				local time_string = subinstr("`time_string'", ":", "_", .)
				local time_string = subinstr("`time_string'", "-", "_", .)
				graph save "rdwa_graph3_`time_string'.gph", replace
				graph export "rdwa_graph3_`time_string'.png", replace width(3000)
				graph export "rdwa_graph3_`time_string'.pdf", replace 
				clear
			}
			restore
		}
		if "`efron2014'"=="on" {
			mata: mata drop `y_star_mata' `t_star_mata' `y_star_t_star_mata' `ones_mata' `y_star_dev_mata' `t_star_mean_mata' `t_star_dev_mata' `z_star_mata' `cov_mata' `bc_mata'
		}
		display " "
		display as text "rdwa finished at $S_TIME  $S_DATE"
		display " "
	}
end

	






