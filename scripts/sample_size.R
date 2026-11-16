# sample_size.R
# simulation of sample size
# June 2023
library(janitor)
library(nimble)
library(broom.mixed)
library(lme4) # for glmer
library(dplyr)
TeachingDemos::char2seed('spain')

# key parameters:
n_policy = 88 - 18 - 15 # journals with a policy, remove journals with no policy date and uncertain date ...
n_no_policy = 176 - 36 - 30 # ... ditto for controls (times 2)
n_journals = n_policy + n_no_policy # total journals
year_start = 2010 # year range for observed data
year_end = 2022
years = year_start:year_end
n_years = year_end - year_start + 1
median_paper_year = 130 - 10 # baseline mean of eligible papers per year per journal (minus 10 for start - see below)
median_paper_year = median_paper_year *2/3 # reduce by one-third due to irrelevant papers
trend_papers = 0 # percent increase in paper numbers over time (kept at zero as can be confounded with policy)
base_rate = 0.5 # base rate (proportion of bad graphs)
base_rate_sd = 0.05 # standard deviation for above (per journal level)
year_slope = 1 # decreasing trend in bad graphs as a percent on absolute scale
efficacy = 1 # percent reduction in bad graphs on an absolute scale due to policy change
n_sim = 1000 # number of simulations
extra_sd = 0.02 # noise for extr-poisson variation

# simulation loop
year.grid = expand.grid(jnum = 1:n_journals, year = year_start:year_end)
sim_res = NULL
for (sim in 1:n_sim){
  
  ## generate the data
  # journals vary in base rate:
  base_per_journal = rnorm(n = n_journals,
                           mean = base_rate, 
                           sd = base_rate_sd )
  # some journals publish more papers than others:
  papers_per_journal = round(10 + rexp(n = n_journals, # use exponential for skew, 10 is minimum
                                       rate = 1/median_paper_year))
  
  # randomly generate policy year
  policy = data.frame(jnum = 1:n_journals,
                      base = base_per_journal) %>%
    mutate(policy_year = round(runif(n=n(), min = 2017, max = 2021)), # most policies came in 2017 to 2018
           policy_year = ifelse(jnum > n_policy, 9999, policy_year)) # dummy year for journals without a policy

  # create simulated data
  frame = data.frame(mean = papers_per_journal) %>%
    mutate(jnum = 1:n()) %>%
    right_join(year.grid, by='jnum') %>% # frame of all journals in all years
    mutate(
      mean = mean *(1+((trend_papers/100))*(year-year_start+1)), # trend in paper numbers over calendar time
      papers = rpois(n = n(), lambda = mean)) %>% # number of papers per journal per year
    left_join(policy, by='jnum') %>%
    mutate(has_policy = as.numeric(year > policy_year), # binary pre-post policy
           yearc = year - 2015, # centre year to help with model convergence
           noise = rnorm(n=n(), mean = 0, sd = extra_sd), # add extra-poisson variation - causing false positive issues
           p = base - (year-2010)*(year_slope/100) - (has_policy*efficacy/100) + noise, # probability of bad graphs
           p = pmax(p, 0),
           observed = rbinom(n = n(), size = papers, prob = p)) # randomly generated bad paper
  
  # binomial version, gives same answer
  # tryCatch for occasional error
  # modelb = tryCatch(glmer(cbind(observed, papers-observed) ~ yearc + has_policy + (1|jnum), 
  #               family = binomial,
  #               data = frame), 
  #               error = function(e) print('Convergence error'))
   modelp = tryCatch(glm(observed ~ yearc + has_policy + offset(log(papers)) + factor(jnum), # fixed effect for journal
                           family = quasipoisson,
                           data = frame), 
                     error = function(e) print('Convergence error'))
   if(class(modelp)[1] == 'character'){next} # skip if error
   ests = tidy(modelp) %>%
     filter(term %in% c('has_policy','yearc'))
  
  ## bayesian model - too slow!
  # source('bayes_model.R')
  # ests = as.data.frame(mcmc) %>%
  #   clean_names() %>%
  #   tibble::rownames_to_column() %>%
  #   filter(rowname == 'effect')
  
  sim_res = bind_rows(sim_res, ests)
  
  if(sim%%50==0)cat('Simulation', sim, '\n')
  
}
hist(sim_res$estimate) # quick histogram

# power
mutate(sim_res, sig = p.value < 0.05) %>%
  group_by(term) %>%
  summarise(n = n(),
            power = sum(sig)/n)

# visualise one simulation
library(ggplot2)
greys = grey(runif(n = n_journals, min=0.01, max=0.8))
jplot = ggplot(data= frame, aes(x=year, y=papers, group=factor(jnum), col=factor(jnum)))+
  scale_color_manual(values = greys)+
  geom_line()+
  ylab('Number of papers')+
  scale_x_continuous(breaks =years)+
  theme_bw()+
  theme(legend.position = 'none')
jpeg('figures/simulation_journal_numbers.jpg', width=5, height = 4.5, units='in', res=400)
print(jplot)
dev.off()
# 
frame = mutate(frame, p_obs = observed / papers)

pplot = ggplot(data= frame, aes(x=year, y=p, group=factor(jnum), col=factor(jnum)))+
  scale_color_manual(values = greys)+
  geom_line()+
  ylab('Proportion of bad graphs')+
  scale_x_continuous(breaks =years)+
  theme_bw()+
  theme(legend.position = 'none', 
        panel.grid.minor = element_blank())+
  facet_wrap(~policy_year)
jpeg('figures/simulation_observed_proportions.jpg', width=5, height = 4.5, units='in', res=400)
print(pplot)
dev.off()
