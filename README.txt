Title: R code for fitting multi-state Markov model and producing figures in 'Quantification of the natural history of visceral leishmaniasis and consequences for control'

Version: 0.1.0

Author: Lloyd Chapman

Email: Lloyd.Chapman@lshtm.ac.uk

Description: R code for fitting 5-state and 6-state Markov models of the natural history of visceral leishmaniasis (VL) described in [1], using the msm R package [2]. The main files for defining and fitting the models are 'run_msm_BNVL2004data_all_deaths_and_relapse.R' and 'run_msm_BNVL2004data_extra_asympt_state_all_deaths_and_relapse.R'. The repository also contains code for plotting Kaplan-Meier survival curves for progression to VL according to seropositivity and seroconversion status ('plot_KM_curve_for_seropstvty2.R' and 'plot_KM_curve_for_serocnvsn2.R'). The code outputs the results shown in Fig. 3, equation (2), Tables 3 and 4, and Additional File 1 of [1].

Repository contents:
createSurvivalFrame.R
hzd_ratio_with_conf_intvls.R
pair_log_rank_test.R
plot_KM_curve_for_serocnvsn2.R
plot_KM_curve_for_seropstvty2.R
plot_obs_num_in_each_state.R
plot_prevalence.R
plot_survival_prob.R
process_BNVL2004data_all_deaths_and_relapse.R
process_BNVL2004data_extra_asympt_state_all_deaths_and_relapse.R
run_msm_BNVL2004data_all_deaths_and_relapse.R
run_msm_BNVL2004data_extra_asympt_state_all_deaths_and_relapse.R
(See individual files for a description of their function.)

Developed in: R 3.2.1. R Foundation for Statistical Computing, Vienna, Austria, 2015. https://www.r-project.org/. 

Installation: R, which is free to install, must be installed to run the code. It can be downloaded from https://www.r-project.org/. Installing RStudio as an R development environment is recommended. RStudio can be downloaded at https://www.rstudio.com/products/rstudio/download/. N.B. The code cannot be run without the associated dataset, which must be requested from the author. The code uses the multi-state modelling R package msm, which can be installed by typing "install.packages("msm")" at the command prompt in RStudio. Clone/download and save the contents of this repository in a folder, open R/RStudio and change the working directory to the folder in which the files and data are saved. The code for fitting the 5-/6-state Markov models can then be run by typing "source('run_msm_BNVL2004data_all_deaths_and_relapse')" or "source('run_msm_BNVL2004data_extra_asympt_state_all_deaths_and_relapse.R')" at the R command prompt. 

License: GNU Affero General Public License v3.0 (http://www.gnu.org/licenses/agpl-3.0.txt)

References: [1] Chapman LAC, et al. Quantification of the natural history of visceral leishmaniasis and consequences for control. Parasites & Vectors, 2015; 8:521. https://doi.org/10.1186/s13071-015-1136-3
[2] Multi-state modelling with R: the msm package. MRC Biostatistics Unit. 2014. https://cran.r-project.org/web/packages/msm/vignettes/msm-manual.pdf
