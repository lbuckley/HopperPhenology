## Zipkin_et_al_Ecology_2014
*Zipkin E.F., Thorson J.T., See K., Lynch H.J., Grant E.H.C., Kanno Y., Chandler R.B., Letcher B.H., and Royle J.A. 2014. Modeling structured population dynamics using data from unmarked individuals. Ecology. 95: 22-29.*

These data and scripts are associated with Zipkin et al. 2014 Ecology.

**Data:**

saldata.aggregated_A.csv - Adult salamander counts for 21 sites from 2005-2011, with two sampling events per year
saldata.aggregated_J.csv - Juvenile salamander counts for 21 sites from 2005-2011, with two sampling events per year

**Code:**

Wrapper R code_salamander application.R - Start with this file. This code reads and formats the data and then fits a structured Dail-Madsen model using the associated JAGS code "JAGS code_salamander application.R". The code produces a summary of the parameter estimates for the model and model diagnostics (e.g., traceplots and R-hat statisic).

salamander_leslie matrix.R" - estimates the population growth rate, lambda, using the full posterior distribution.
