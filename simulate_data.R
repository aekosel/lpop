######################################################################
###Author: aekosel
###Purpose: Because the dataset used in the paper is not publicly
###available, we simulate a similar dataset for the reader to use.
######################################################################

#create discrete uniform generator
rdunif <- function(n, min, max){
	sample(min:max, n, replace=T)
}

#set number of subjects
n = 5000

#create example variables
set.seed(123456)
age = c(rdunif(floor(2*n/3), 70, 85), rdunif(ceiling(n/3), 70, 97))
disability_score_base = trunc( rnorm(n) * 6 + 12 )
disability_score_base = ifelse(disability_score_base < 0, 0, 
                               disability_score_base)
disability_score_base = ifelse(disability_score_base > 24, 24, 
                                disability_score_base)
#
recovered = rbinom(n, 1, .3)
disability_score_1yr = recovered*trunc( (runif(n)^2)*5 )  +
	(1 - recovered)*( trunc(  2 + 
                    0.85 * disability_score_base +
                    0.20 * (age-75) + rnorm(n)*5 ) )
disability_score_1yr = ifelse(disability_score_1yr < 0, 0, 
                               disability_score_1yr)
disability_score_1yr = ifelse(disability_score_1yr > 24, 24, 
                                disability_score_1yr)
#
pain_back_base = trunc( 3 + (6/24) * disability_score_base + rnorm(n)*2 )
pain_back_base = ifelse( pain_back_base < 0, 0, 
                               pain_back_base )
pain_back_base = ifelse( pain_back_base > 10, 10, 
                               pain_back_base )
#
pain_leg_base = trunc( 4 + (3/24) * disability_score_base + rnorm(n)*2 )
pain_leg_base = ifelse( pain_leg_base < 0, 0, 
                               pain_leg_base )
pain_leg_base = ifelse( pain_leg_base > 10, 10, 
                               pain_leg_base )
#
pred_score = fitted( lm( disability_score_base ~ pain_back_base +
                                                 pain_leg_base ) )
#
#save dataset
sim_data = cbind(age,
	disability_score_base, disability_score_1yr,
	pain_back_base, pain_leg_base,
	pred_score)

write.csv(sim_data, "simulated_data.csv", row.names=F)
