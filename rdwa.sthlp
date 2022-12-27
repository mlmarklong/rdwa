{smcl}
{* *! version 1.0.1 27dec2022}{...}
{title:Title}

{phang}
{bf:rdwa} {hline 2} Regression(s) discontinuity weighted average estimator

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:rdwa}
varlist(min=2
max=2)
[{help if}]
[{help in}]
[{cmd:,}
{it:options}]

{pstd}
where {it:varlist} is 
{p_end}
		{it:Y} {it:X}

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}

{synopt:{opt c(#)}}  		Cutoff value for receipt of treatment. Default = 0.{p_end}
{synopt:{opt p_min(#)}}  	# = Minimum polynomial order to be considered. Default = 0.{p_end}
{synopt:{opt p_max(#)}}  	# = Maximum polynomial order to be considered. Default = 4.{p_end}
{synopt:{opt rbc}}  		Generate bias-corrected RDTE estimates.{p_end}
{synopt:{opt kernel(string)}}  	Kernel used for weighting regressions ("uni" (default) or "tri").{p_end}
{synopt:{opt samples(#)}}  	# of bootstrapped samples used to generate regressions. Default = 200.{p_end}
{synopt:{opt efron2014}}  	Include Efron's (2014) method for standard error estimates.{p_end}
{synopt:{opt details_off}} 	Do no show estimates from each bootstrapped sample.{p_end}
{synopt:{opt graph1}}  		Generate and save a figure showing the distribution of estimated RDTEs.{p_end}
{synopt:{opt graph2}}  		Generate and save a figure showing the effective weight on obs.{p_end}
{synopt:{opt graph3}}  		Generate and save a scatterplot with regression lines.{p_end}
{synopt:{opt graph3_lines(#)}} 	# of regression lines shown in graph 3. Default = min(200, samples).{p_end}

{synoptline}
{p2colreset}{...}
{p 4 6 2}

{marker description}{...}
{title:Description}

{pstd}
{cmd:rdwa} does the following: 
(a) presents estimates of the regression discontinuity treatment effect (RDTE) assuming a local linear specification with 
bandwidths determined by the methods of Calonico, Cattaneo, and Titiunik (CCT, 2014a);
(b) selects the polynomial order using the method of Pei, Lee, Card, and Weber (PLCW, 2021) with bandwidths 
determined by the methods of Calonico, Cattaneo, and Titiunik (CCT, 2014a); and
(c) uses the method of Long and Rooklyn (2022), which generates 200 (default) estimates of the RDTE 
using 200 bootstrapped samples and applying the PLCW algorithm to each sample, and presents the average of the 
RDTE estimates and an estimate for the standard deviation of the average of the 
RDTE estimates. For more information on this algorithm and its wisdom and efficacy, see 
{browse "http://evans.uw.edu/profile/long":{it:Regression(s) Discontinuity}.}

{pstd} Using the specifications selected by the Long and Rooklyn algorithm, the {bf:rdwa} command then generates a 
spiffy graph demonstrating both the PLCW chosen "best" specification and the variance in the estimated RDTE that 
comes from uncertainty about what is the best polynomial order and the best bandwidth.

{pstd} {bf:User Notes:}

{pstd} (0) The Long and Rooklyn algorithm bootstraps the PLCW method, which builds on the CCT method. To implement
{bf:rdwa}, it is necessary for the user to obtain the corresponding PLCW and CCT programs using the following commands:
{stata `"net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace"' : net install rdmse, from(https://raw.githubusercontent.com/peizhuan/rdmse/master) replace}
and {stata `"ssc install rdrobust, replace"' : ssc install rdrobust, replace}

{pstd} (1) The first variable in the varlist is the outcome of interest (Y). The second variable in the varlist is the 
running variable (X). No other variables should be included. (At this time, the {bf:rdwa} command is not set-up to 
handle covariates).

{pstd} (2) "c(#)" is the cutoff for receiving treatment. This program produces an RD estimate assuming a strict cutoff 
in treatment at the threshold.  (At this time, the {bf:rdwa} command is not set-up to estimate fuzzy RDs or RD kink designs).
It is assumed that an observation with X=c is on the right side of the discontinuity. If an observation that is exactly
at the cutoff should be considered to be with the group on the left side of the discontinuity, then set the 
cutoff to be cutoff+epsilon, where epsilon is a very small number that will successfully cleave the treatment and control 
groups into two appropriate parts. (Alternatively, generate a new variable that is -X (i.e., reverse the order of X) and 
apply the {bf:rdwa} command using this transformed variable).

{pstd} (3) "p_min(#)" ("p_max(#)") defines the minimum (maximum) order of the polynomial that you wish to be considered as a 
candidate specification. p_min has 0 as its default, i.e., a kernel-weighted regression of y on "treated" for observations 
within the bandwidth, where "treated"=1 if X>=cutoff. p_max has 4 as its default, i.e., a kernel-weighted regression 
of Y on X, X^2, X^3, X^4, treated, treated*X, treated*X^2, treated*X^3, and treated*X^4 for observations within the bandwidth. 

{pstd} (4) The default estimates are "conventional" and use conventional standard errors. If the user wishes to generate
bias-corrected RDTE estimates, using the methods of CCT (2014a), for both the PLCW and Long and Rooklyn procedures then the 
user should use the "rbc" option. This option will produce robust standard errors for the PLCW results, again using the methods 
of CCT (2014a). For the Long and Rooklyn results, the standard error is estimated using the procedure outlined in Efron (2014) 
from Equation 3.6. 

{pstd} (5) The default kernel weight is uniform (i.e., equal weight for observations within the bandwidth rather than 
triangular (i.e., with more weight placed on observations near to the threshold and with a steady decline throughout the range 
of X on which the regression is run). As documented in our associated paper, a uniform weighting scheme might be better for our 
algorithm as the average RDTE coming from various bandwidths provides less weight on observations away from the cutoff value, 
even when using uniform weights for each individual regression. If the RDTE estimates were based on a single regression, as 
standard in the RD literature, then triangular kernels would be recommended per the findings of Cheng, Fan, and Marron (1997), 
who conclude that for local polynomial estimation the triangular weight function is best for any particular order. 
For the interested user, {bf:rdwa} allows triangular weights to be used ("kernel(tri)"). 

{pstd} (6) "samples(#)" is the number of bootstrapped samples that are drawn to generate a set of estimated RDTEs which are 
then averaged to generate an estimate of the RDTE. The default value is 200 and thus the algorithm will generate 200 
estimates. If the user selects a higher value for "samples", then more than 200 estimates will be generated. A larger value of 
samples means slower time to completion of the command, but more accurate estimates of the RDTE. The user might want to
set samples to a low number as a preliminary run to see what the algorithm produces (e.g., "10" in the example below).

{pstd} (7) "efron2014" generates the bias-corrected standard deviation for the smoothed bootstrap estimate (i.e., Efron's 
Equation 7.25 on p. 1006). Efron shows that in the "ideal case" (i.e., when the number of bootstrap samples equals N^N and with 
each possible combination of the N observations chosen once), the standard deviation for the smoothed bootstrap estimate 
(i.e., Equation 3.4 on p. 995) is less than the standard deviation of the various bootstrap estimates (i.e., Equation 2.4 on p. 993). 
That is, the "standard confidence interval" that is constructed using the standard deviation of the various bootstrap estimates
(which is the default method used in this program) is too wide. Efron's method corrects this issue. Note, however, that Efron's method
for estimating the standard deviation in the "ideal case" still requires a high number of bootstrapped samples.  Efron uses 4,000, 
which is far more than the 200 used as the default setting in this program. When Efron's method is applied using a small number of 
bootstrapped samples, the resulting estimate of the standard deviation for the smoothed bootstrap estimate is much noisier than
that found using the standard deviation of the various bootstrap estimates and it may incorrectly suggest that the "standard 
confidence interval" is too narrow. The user is thus cautioned that this option is best used with a high setting of the samples option, 
e.g., samples(4000).

{pstd} (8) {bf:rdwa} shows a summary of the results from the 200 bootstrapped samples. If you do want to see the individual
estimates of the RDTE for each bootstrapped sample, then you should use the "details_off" option.

{pstd} (9) The "graph1", "graph2", and "graph3" options generate and save three graphs: rdwa_graph1_####_##_##_##_##_##.&&&, 
rdwa_graph2_####_##_##_##_##_##.&&&, and rdwa_graph3_####_##_##_##_##_##.&&&, where ####_##_##_##_##_## is a record of the
year_month_day_hour_minute_second when the graph was saved and &&& is gph, png, and pdf. "graph1" shows the distribution of the 
estimated RDTEs from the 200 bootstrapped samples. "graph2" shows the effective weight placed on observations to generate the 
RDTEs (with 100% denoting that observations with this value of X are inside the bandwidth in each of the 200 specifications). 
"graph1" and "graph2" are not produced if the "samples" option is set equal to 1. "graph3"
shows the raw scatterplot of the data and the 200 regression lines that correspond to the 200 specifications. Graphs are saved
in three forms: (a) Stata's ".gph" format, which can be opened and manipulated (e.g., colors changed) using Stata's "Graph Editor"; 
(b) a moderately high resolution ".png" format, which can be easily pasted into a Word or other document; and (c) a high
resolution ".pdf" format.  Other formats can be created by opening the ".gph" and using "graph export...". Of these, graph3 is 
likely to be the most useful for most user.  graph3 relies on the latest version of Ben Jann's 
{bf:addplot} command. This command must be installed prior to producing the {bf:rdwa} graph using the following command:
{stata `"net install addplot, replace from(https://raw.githubusercontent.com/benjann/addplot/main/)"' : net install addplot, replace from(https://raw.githubusercontent.com/benjann/addplot/main/)}

{pstd} (10) The default setting for the "graph3_lines" option = min(200, "samples"). If using the default setting for "samples" (200), 
the graph will show all 200 regression lines from these 200 bootsrapped samples. If the user sets a value for "samples" that 
is greater than the user-set (or default) value of "lines", then the lines shown will be taken from the first "lines" bootstrapped
samples. For example, if the user includes the following options, "samples(500) lines(100)", then the algorithm will take 500
bootstrapped samples, generate estimates of the RDTE for each using the methods of PLCW(2021), average these 500 estimates to 
generate the overall estimate of the RDTE, and then generate a graph that shows the first 100 out of 500 regression lines. Note
that the graph adjusts the extent of transparency of each drawn regression line. With the default value of 200 for "lines",
each line has a 5% level of color -- if 20 regression lines overlap at a given point, the full color will be shown (as 5%*20
= 100%).  

{pstd} (11) 
{title:Example}

{pstd}For reproducable results, clear all and set seeds. {p_end}

{col 9}{stata `"clear all"' : clear all}
{col 9}{stata `"set seed 123456789"' : set seed 123456789}

{pstd}Imagine a fictional treatment applied to indviduals with a body mass index (i.e., weight in kilograms / height in meters squared) >= 30.
Suppose that this treatment is intended to lower serum cholesterol (mg/dL). Does this treatment have an effect at the 30-BMI threshold? (While the 
point estimates of the RDTE for this fictional treatment suggest a negative effect on serum cholesterol and this pseudo-effect, it is not found to be 
statistically significant). {p_end}

{col 9}{stata `"use tcresult height weight in 1/500 using http://www.stata-press.com/data/r15/nhanes2.dta"' : use tcresult height weight in 1/500 using http://www.stata-press.com/data/r15/nhanes2.dta}
{col 9}{stata `"gen bmi=weight/(height/100)^2"' : gen bmi=weight/(height/100)^2}
{col 9}{stata `"label var bmi "Body Mass Index""' : label var bmi "Body Mass Index"}
{col 9}{stata `"rdwa tcresult bmi, c(30) efron2014 graph1 graph2 graph3 samples(50)"' : rdwa tcresult bmi, c(30) efron2014 graph1 graph2 graph3 samples(50)}

{title:Stored results}
{synoptset 15 tabbed}{...}
{p2col 5 15 19 2: Scalars}{p_end}
{synopt:{cmd:e(FULL_obs_l)}}  {p_end}
{synopt:{cmd:e(FULL_obs_r)}}  {p_end}
{synopt:{cmd:e(CCT_P)}}  {p_end}
{synopt:{cmd:e(CCT_BW)}}  {p_end}
{synopt:{cmd:e(CCT_RDTE)}}  {p_end}
{synopt:{cmd:e(CCT_RDTE_se)}}  {p_end}
{synopt:{cmd:e(CCT_95ci_lower)}}  {p_end}
{synopt:{cmd:e(CCT_95ci_upper)}}  {p_end}
{synopt:{cmd:e(CCT_obs_l)}}  {p_end}
{synopt:{cmd:e(CCT_obs_r)}}  {p_end}
{synopt:{cmd:e(PLCW_P)}}  {p_end}
{synopt:{cmd:e(PLCW_BW)}}  {p_end}
{synopt:{cmd:e(PLCW_RDTE)}}  {p_end}
{synopt:{cmd:e(PLCW_RDTE_se)}}  {p_end}
{synopt:{cmd:e(PLCW_95ci_lower)}}  {p_end}
{synopt:{cmd:e(PLCW_95ci_upper)}}  {p_end}
{synopt:{cmd:e(PLCW_obs_l)}}  {p_end}
{synopt:{cmd:e(PLCW_obs_r)}}  {p_end}
{synopt:{cmd:e(LR_P0_freq)}}  {p_end}
{synopt:{cmd:e(LR_P1_freq)}}  {p_end}
{synopt:{cmd:e(LR_P2_freq)}}  {p_end}
{synopt:{cmd:e(LR_P3_freq)}}  {p_end}
{synopt:{cmd:e(LR_P4_freq)}}  {p_end}
{synopt:{cmd:e(LR_ave_BWs)}}  {p_end}
{synopt:{cmd:e(LR_sd_BWs)}}  {p_end}
{synopt:{cmd:e(LR_min_BWs)}}  {p_end}
{synopt:{cmd:e(LR_max_BWs)}}  {p_end}
{synopt:{cmd:e(LR_ave_obs_l)}}  {p_end}
{synopt:{cmd:e(LR_ave_obs_r)}}  {p_end}
{synopt:{cmd:e(LR_min_obs_l)}}  {p_end}
{synopt:{cmd:e(LR_min_obs_r)}}  {p_end}
{synopt:{cmd:e(LR_max_obs_l)}}  {p_end}
{synopt:{cmd:e(LR_max_obs_r)}}  {p_end}
{synopt:{cmd:e(LR_ave_RDTEs)}}  {p_end}
{synopt:{cmd:e(LR_pct_95ci_lower)}}  {p_end}
{synopt:{cmd:e(LR_pct_95ci_upper)}}  {p_end}
{synopt:{cmd:e(LR_pct_p_1tail)}}  {p_end}
{synopt:{cmd:e(LR_pct_p_2tail)}}  {p_end}
{synopt:{cmd:e(LR_sd_ave_RDTEs)}}  {p_end}
{synopt:{cmd:e(LR_normal_p_2tail)}}  {p_end}
{synopt:{cmd:e(LR_normal_95ci_lower)}}  {p_end}
{synopt:{cmd:e(LR_normal_95ci_upper)}}  {p_end}
{synopt:{cmd:e(LR_efron2014_sd_ave_RDTEs)}}  {p_end}
{synopt:{cmd:e(LR_efron2014_normal_p_2tail)}}  {p_end}
{synopt:{cmd:e(LR_efron2014_normal_95ci_lower)}}  {p_end}
{synopt:{cmd:e(LR_efron2014_normal_95ci_upper)}}  {p_end}

{title:Acknowledgments}
{p}

{pstd} Support for this research came from the U.S. Department of Education’s Institute of Education Sciences (R305A140380) and a 
Eunice Kennedy Shriver National Institute of Child Health and Human Development research infrastructure grant (R24 HD042828) to the 
Center for Studies in Demography & Ecology at the University of Washington. Helpful comments were provided by Tamre Cardoso, Wenyu 
Chen, Brian Dillon, Ariane Ducellier, Dan Goldhaber, Aureo de Paula, Laura Peck, Jon Smith, Sarah Teichman, Seth Temple, Jake Vigdor, 
Ted Westling, Xiaoyang Ye, and Association for Public Policy and Management conference and University of Washington seminar audience 
members. Excellent research assistance was provided by Ben Glasner and Tom Lindman. Finally, we thank 
Ben Jann for helpful comments regarding Ben's "addplot" command.

{title:Citation}
{p}

{pstd} Thank you for citing {bf:rdwa} as follows: Mark C. Long and Jordan Rooklyn. 2022. 
"rdwa: Stata command for the algorithm described in 'Regression(s) Discontinuity'".

{title:Author}
{p}

{phang}
Mark C. Long, University of Washington (corresponding author). Email {browse "mailto:marklong@uw.edu":marklong@uw.edu}{p_end}
{p2colreset}{...}
{phang}
Jordan Rooklyn, City of Talent{p_end}

{title:References}

{phang}
Long, M., and J. Rooklyn. 2022.
{browse "http://evans.uw.edu/profile/long":{it:Regression(s) Discontinuity: Using Bootstrap Aggregation to Yield Estimates of RD Treatment Effects}}.
Working paper.
{p_end}

{phang}
Pei, Z., D.S. Lee, D. Card, and A. Weber. 2021. Local Polynomial Order in Regression Discontinuity Designs. 
{browse "https://doi.org/10.1080/07350015.2021.1920961":{it:Journal of Business & Economic Statistics}}
{p_end}

{phang}
Calonico, S., M.D. Cattaneo, and R. Titiunik. 2014a. Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs.
{browse "https://www.jstor.org/stable/43616914":{it:Econometrica}}, 82(6): 2295-2326.

{phang}
Calonico, S., M.D. Cattaneo, and R. Titiunik. 2014b. Robust Data-Driven Inference in the Regression-Discontinuity Design
{browse "https://journals.sagepub.com/doi/pdf/10.1177/1536867X1401400413":{it:Stata Journal}}, 14(4): 909-946.

{phang}
Efron, B. 2014. Estimation and Accuracy After Model Selection.
{browse "https://doi.org/10.1080/01621459.2013.823775":{it:Journal of the American Statistical Association}}, 
109(507):991-1007.
{p_end}

{phang}
Cheng, M., Fan, J, Marron, J.S. 1997. On Automatic Boundary Corrections.
{browse "https://www.jstor.org/stable/2959068":{it:Annals of Statistics}}
, 25(4):1691-1708. 
{p_end}

{title:Disclaimer and Warning}

{pstd} This software is provided "as is" without warranty of any kind, either 
expressed or implied, including, but not limited to, the implied warranties 
of merchantability and fitness for a particular purpose. The entire risk as 
to the quality and performance of the program is with you. Should the 
program prove defective, you assume the cost of all necessary servicing, 
repair, or correction. In no event will the copyright holders or their 
employers, or any other party who may modify and/or redistribute this 
software, be liable to you for damages, including any general, special, 
incidental, or consequential damages arising out of the use of this software.  
Ability to use the program (including but not limited to loss of data or 
data being rendered inaccurate or losses sustained by you or third parties or
a failure of the program to operate with any other programs), even if such 
holder or other parties are advised of the possibility of such damages.




