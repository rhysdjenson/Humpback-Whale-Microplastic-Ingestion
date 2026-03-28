library (TruncatedNormal)
library(boot)

#### This model converts the number of MP per gram of feces (dry) to the number of MP per gram prey (wet). 
#### This value is multiplied by estimated mass of food ingested daily.
#### For whale eating 30% juvenile herring and 70% krill **** Based on WGS results from Ch.4

set.seed(100)
# simulate dry weight percent of krill and herring
herringdw <- rnorm(1000, mean = 20, sd = 2.4)
krilldw <- rnorm(10000, mean = 18.97, sd = 1.42)
# weighted sampling (70% krill, 30% herring)
combined <- c(
  sample(krilldw, size = 70, replace = TRUE),
  sample(herringdw, size = 30, replace = TRUE)
)
mean.dwprey <- mean(combined)
sd.dwprey   <- sd(combined)

# Manganese values (PPM) for freeze-dried feces and prey.
FecalMn_Avg <- 32.54
FecalMn_SD <- 2.79894996
KrillMn_Avg <- 5.19
KrillMn_SD <- 1.466162213
HerringMn_Avg <- 1.61

#Weighted average (70% krill, 30% herring) for manganese values between both prey types: 
#Total mean
totalpreyMn_Avg <- 0.7 * KrillMn_Avg + 0.3 * HerringMn_Avg
# total variance
totalpreyMn_Var <- 0.7*(KrillMn_SD^2 + (KrillMn_Avg - totalpreyMn_Avg)^2) +
  0.3*(0 + (HerringMn_Avg - totalpreyMn_Avg)^2)
#Total SD
totalpreyMn_SD <- sqrt(totalpreyMn_Var)




# Mean and SD for the ratio of feces to prey (d.w.)
ratio_with_error <- function(mean_x, sd_x, mean_y, sd_y) {
  mean_z <- mean_x / mean_y
  rel_error <- sqrt(
    (sd_x / mean_x)^2 +
      (sd_y / mean_y)^2
  )
  sd_z <- mean_z * rel_error
  
  list(mean = mean_z, sd = sd_z)
}
result<-ratio_with_error(totalpreyMn_Avg, totalpreyMn_SD,FecalMn_Avg, FecalMn_SD)
dryratio_fecesprey_avg<-result$mean 
dryratio_fecesprey_sd<-result$sd 

# Conversion ratio (equation 3.3 and 3.4 from Ch.3)
fecalMPdry <- rtnorm(1000, mu=24.36, sd=16.37, lb=10, ub=68.75)
PreyFecesRatio <- rtnorm(1000, mu=dryratio_fecesprey_avg, sd=dryratio_fecesprey_sd, lb=0, ub=Inf)
MPpergrampreywet<-(fecalMPdry*PreyFecesRatio)*(rnorm(1000, mean = mean.dwprey, sd = sd.dwprey)/100)
mean(MPpergrampreywet)

# Next, estimate consumption. Formula From Witteveen et al., 2006
m<-30000
K <- 0.88*0.3 + 0.74*0.7 
E<- 192*(m^0.75)
I <- (E/K)*(1/1000)
r <- I*1000

((r/1000)/m)*100 # Whales feeding on 70% euphausids and 30% herring eat 1.87% bw per day


MPingestedperday<-r*MPpergrampreywet
mean(MPingestedperday)
median(MPingestedperday)
hist(MPingestedperday)
quantile.95MPingestedperday<-quantile(MPingestedperday, probs = c(0.05, 0.95))


#Bootstrap resampling (n = 1000) yielded a mean ingestion estimate of 418,227 MPs/day 
#with low bias and a standard error of 12,298 MPs/day. The 95% percentile confidence interval 
#ranged from 396,731 to 424,217 MPs/day, indicating relatively low uncertainty and high stability of 
#the estimate.

boot_mean <- function(data, i) mean(data[i], na.rm = TRUE)
boot_res <- boot(MPingestedperday, boot_mean, R = 1000)
bootresult<-boot.ci(boot_res, type = "perc")
conf<-bootresult$percent
ci_lower <- conf[1, ncol(conf)-1]
ci_upper <- conf[1, ncol(conf)]

# 30,000 kg HB feeding on %70 krill and 30% juv herring  may ingest mean 418227.2 MP a day (95% CI, 401190.5-436171.9).
#Mass ingested estimate
# Density of cellulose = 1.5 g/cm³ = 1500 kg/m³
# Average MP isolated was 900  micron x 30 micron, 72 % cellulose.
# Average 0.95 micrograms per MP
#Upper quantile of highest scanrio (100% krill) = 1,094,463.5 / day

meanmassingested<-(9.54*10^-7)*mean(MPingestedperday) # mean = 0.3989887 g/day
(9.54*10^-7)*ci_lower #5% = 0.3827357  g/day
(9.54*10^-7)*ci_upper #95% =  0.416108 g/day

# ~ mean 0.3989887 grams of MPs per day (95CI:0.382735-0.416108)

((meanmassingested/1000)/30000)*100

# Represents 1.329962e-06 % of 30,000kg mass consumed daily... Insignificant.
