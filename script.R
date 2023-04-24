library(tidyverse)
library(here)
library(rio)
library(janitor)
library(brms)
library(tidybayes)

seperate()
septic_shock_90 <- import_list(here("data/sepsis_data.xlsx"))$`90-day septic shock` |>
  clean_names() |>
  pivot_longer(cols = contains("events"),names_to = "group",values_to = "n") |>
  separate_wider_delim(n,delim = "/", names = c("events","total")) |>
  mutate(across(group,~str_remove(.x,"_.*")),
         year = str_extract(study, "\\d{4}"),
         across(study,as.factor),
         across(contains(c("events","total","year")),as.numeric))

model1 <- brm(events | trials(total) ~ 0 + Intercept + group + (1 + group || study), # multilevel (AKA random effects) logistic model, allowing intercept (base event rate) and difference to vary between studies (|| rather than | excludes covariance between the two, which can't be estimated here)
               family = binomial(link ="logit"), # binomial likelihood, with logit link. In brms, these models are specified using "nubmer of events" | trials("number of subjects/trials") on the left hand side of the formula
               prior = prior(normal(0,1.5),class = "b",coef="Intercept") + #prior for base rate, different than original analysis (more sensible)
                 prior(normal(0,2.82),class = "b") + # prior for difference propofol - comparator on log-odds scale, same as original analysis
                 prior(cauchy(0,0.5), class = "sd", group = "study", coef = "Intercept") +
                 prior(cauchy(0,0.2), class ="sd"), # prior for standard deviation for random effects, same as original analysis
               data = septic_shock_90, #original data
               control = list( adapt_delta = 0.95),
               file = here("fits/model1.rds"), #saves the model fit
               file_refit ="on_change")

or_data <- model1 |>
  spread_rvars(b_groupcrystalloid) |>
  transmute(or=exp(-b_groupcrystalloid)) |>
  mutate(model="")

#comparison labels
or_labels <-
  or_data |>
  mutate(pp = mean(or < 1)) |> #compute posterior probability for harm from propofol by simply calculating what proportion of the draws have an rr above 1
  median_hdci(or) |> # takes .epred (the rr) and computes the median and 95% highest density interval (credible interval)
  mutate(across(where(is.numeric),  ~ round(.x, 2)))

# plot
or_data|>
  ggplot()+
  aes(dist=or,
      fill=after_stat(x>1),
      y=model) + #after_stat to give the fill colour difference
  stat_halfeye()+ # nice geom for baysian samples, density + interval plot
  geom_text(aes(label = glue::glue("{or} [{.lower}, {.upper}]\n Post. prob. {pp}"),x=Inf),data=or_labels,hjust = "inward") +
  geom_vline(xintercept = 1,linetype = 2) +
  scale_x_continuous(name = "Odds ratio") +
  scale_y_discrete(name = "") +
  scale_fill_manual( values = c("steelblue","pink")) +
  theme_linedraw() +
  theme(legend.position = "none") + #remove legend
  labs(caption = "@load_dependent")

