---
title: "R Notebook"
output:
  html_document:
    keep_md: true

---
A model to calculate the optimal number of workers in an 8-hour shift. The workers change shifts either at the beginning of the shift or at the middle, so that at every given shift workers from the current shift and the previous shift are present.
We know the mean value of the required number of workers, but this requirement is stochastic. We suppose that it follows a Poisson distribution. For example if a certain 4 hour window requires 9 workers, then 


```r
cat('95% of the time the required workers will be less than', qpois(0.95,9),
    'and more than', qpois(0.95,9,lower.tail = F))
```

```
## 95% of the time the required workers will be less than 14 and more than 4
```
Assume that we have the following mean rates and we set the constraint so that 60% of the time they will be satisfied.

```r
mean_rates = c(9,22,25,17,30,13)
CHANCECTR = 0.6 #60% of the times the requirements will be satisfied. The higher the bigger the value.
```




```r
library(data.table)
library(lpSolveAPI)
library(kableExtra)
library(lubridate)
```

```
## 
## Attaching package: 'lubridate'
```

```
## The following objects are masked from 'package:data.table':
## 
##     hour, isoweek, mday, minute, month, quarter, second, wday, week,
##     yday, year
```

```
## The following objects are masked from 'package:base':
## 
##     date, intersect, setdiff, union
```

```r
library(stringr)
options(warn=-1)

DF = data.frame(val = 1:6,
                dt_start = paste0(as.character(seq(0,20,by=4)),':','00')  ,
                dt_end = paste0(as.character(seq(4,24,by=4)),':','00') ,
                effective_rate_ctr = sapply(mean_rates, function(x) qpois(CHANCECTR,x)),
                mean_rates = mean_rates,
                conf_95_lower = sapply(mean_rates, function(x) qpois(0.95,x)),
                conf_95_uper = sapply(mean_rates, function(x) qpois(0.95,x, lower.tail=F)),
                stringsAsFactors = FALSE)

DF %>%
  kbl() %>%
  kable_styling() %>%
  column_spec(4, background = "red")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> val </th>
   <th style="text-align:left;"> dt_start </th>
   <th style="text-align:left;"> dt_end </th>
   <th style="text-align:right;"> effective_rate_ctr </th>
   <th style="text-align:right;"> mean_rates </th>
   <th style="text-align:right;"> conf_95_lower </th>
   <th style="text-align:right;"> conf_95_uper </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> 0:00 </td>
   <td style="text-align:left;"> 4:00 </td>
   <td style="text-align:right;background-color: red !important;"> 10 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 4 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> 4:00 </td>
   <td style="text-align:left;"> 8:00 </td>
   <td style="text-align:right;background-color: red !important;"> 23 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 15 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> 8:00 </td>
   <td style="text-align:left;"> 12:00 </td>
   <td style="text-align:right;background-color: red !important;"> 26 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 17 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> 12:00 </td>
   <td style="text-align:left;"> 16:00 </td>
   <td style="text-align:right;background-color: red !important;"> 18 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 11 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:left;"> 16:00 </td>
   <td style="text-align:left;"> 20:00 </td>
   <td style="text-align:right;background-color: red !important;"> 31 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 21 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:left;"> 20:00 </td>
   <td style="text-align:left;"> 24:00 </td>
   <td style="text-align:right;background-color: red !important;"> 14 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 7 </td>
  </tr>
</tbody>
</table>
The variables
$X[i] >= 0, i=1:6$ the number of workers

The constraints: The available staff should be more than the effective constraints
$$AvailStaff(X,i) >= VaR(RequiredStaff;ChanceCtr),  i=1:6$$

The objective
$$min(sum_i X[i])$$





```r
initialise_model = function(n){
  nVar = n
  nCtr = nVar
  lprec <- make.lp(nCtr,nVar)
  return(list('nVar' = nVar,
              'nCtr' = nCtr,
              'model' = lprec))
}

#Build variables
build_variables = function(nVar){
  VarX_col_idx = 1:nVar
  obj = rep(1,nVar)
  return(list('variables_index' = VarX_col_idx,
              'objecive_coefficient' = obj))
}

#At each 4-hour the workers of the current and previous shift are present
AvailStaff = function(X,i,n){
  ifelse(i>1,X[i] + X[i-1],X[1] + X[n])
}

#Set constraints
build_constraints = function(nCtr, nVar,model){
  row_lb = matrix(NA, nrow = nCtr)
  row_ub = matrix(NA, nrow = nCtr)
  for(i in 1:nVar){
    row_ub[i] = -DF[i,"effective_rate_ctr"]
    if(i>1){
      set.row(model, i,c(-1,-1),
              indices = c(i,i-1))
    }else{
      set.row(model, i,c(-1,-1),
              indices = c(i,nVar))
    }
    
  }
  return(list("row_lb" = row_lb,
              "row_ub" = row_ub))
}


set_controls = function(lpreq, row_lb, row_ub,
                        nCtr,obj, varIdx){
  #Controls
  #set lhs and rhs
  set.constr.value(lpreq,
                   lhs = row_lb,
                   rhs = row_ub, 1:nCtr)
  set.objfn(lpreq, obj)
  set.type(lpreq,
           columns = varIdx, type = c("integer"))
  
}


#solve the model
solve_program = function(lprec){
  solStatusCode =  solve.lpExtPtr(lprec)
  sols <- get.variables(lprec)
  minimizer <- get.objective(lprec)
  
  return(list('status' = solStatusCode,
              'soloutions' = sols,
              'minimiser' = minimizer)
         )
}
```
Solving the program

```r
init_mod = initialise_model(nrow(DF))
variables_objective = build_variables(init_mod$nVar)
constraints = build_constraints(init_mod$nCtr,
                                init_mod$nVar, init_mod$model)

set_controls(init_mod$model, constraints$row_lb,
             constraints$row_ub, init_mod$nCtr,
             variables_objective$objecive_coefficient,
             variables_objective$variables_index)

# write model to model.lp
#write.lp(lprec, "MaintOpt.lp", type = c("lp"))

solution = solve_program(init_mod$model)
```

Calculate the available staff (from the two shifts) and the chances that the day will require no more workers than the solution


```r
available_staff = sapply(1:length(mean_rates),
                         function(x) AvailStaff(solution$soloutions, x, length(mean_rates)) )

chances = sapply(1:length(mean_rates), function(x) ppois(available_staff[x],
                                           mean_rates[x]))

DF$solution = solution$soloutions
DF$available_staff= available_staff
DF$chances = chances


DF %>%
  kbl() %>%
  kable_styling() %>%
  column_spec(4, background = "red") %>%
  column_spec(8, background = "green") %>%
  column_spec(9, background = "green") %>%
  column_spec(10, background = "green")
```

<table class="table" style="margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;"> val </th>
   <th style="text-align:left;"> dt_start </th>
   <th style="text-align:left;"> dt_end </th>
   <th style="text-align:right;"> effective_rate_ctr </th>
   <th style="text-align:right;"> mean_rates </th>
   <th style="text-align:right;"> conf_95_lower </th>
   <th style="text-align:right;"> conf_95_uper </th>
   <th style="text-align:right;"> solution </th>
   <th style="text-align:right;"> available_staff </th>
   <th style="text-align:right;"> chances </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 1 </td>
   <td style="text-align:left;"> 0:00 </td>
   <td style="text-align:left;"> 4:00 </td>
   <td style="text-align:right;background-color: red !important;"> 10 </td>
   <td style="text-align:right;"> 9 </td>
   <td style="text-align:right;"> 14 </td>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:right;background-color: green !important;"> 10 </td>
   <td style="text-align:right;background-color: green !important;"> 10 </td>
   <td style="text-align:right;background-color: green !important;"> 0.7059883 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 2 </td>
   <td style="text-align:left;"> 4:00 </td>
   <td style="text-align:left;"> 8:00 </td>
   <td style="text-align:right;background-color: red !important;"> 23 </td>
   <td style="text-align:right;"> 22 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 15 </td>
   <td style="text-align:right;background-color: green !important;"> 25 </td>
   <td style="text-align:right;background-color: green !important;"> 35 </td>
   <td style="text-align:right;background-color: green !important;"> 0.9962450 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 3 </td>
   <td style="text-align:left;"> 8:00 </td>
   <td style="text-align:left;"> 12:00 </td>
   <td style="text-align:right;background-color: red !important;"> 26 </td>
   <td style="text-align:right;"> 25 </td>
   <td style="text-align:right;"> 33 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;background-color: green !important;"> 1 </td>
   <td style="text-align:right;background-color: green !important;"> 26 </td>
   <td style="text-align:right;background-color: green !important;"> 0.6293858 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 4 </td>
   <td style="text-align:left;"> 12:00 </td>
   <td style="text-align:left;"> 16:00 </td>
   <td style="text-align:right;background-color: red !important;"> 18 </td>
   <td style="text-align:right;"> 17 </td>
   <td style="text-align:right;"> 24 </td>
   <td style="text-align:right;"> 11 </td>
   <td style="text-align:right;background-color: green !important;"> 17 </td>
   <td style="text-align:right;background-color: green !important;"> 18 </td>
   <td style="text-align:right;background-color: green !important;"> 0.6549584 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 5 </td>
   <td style="text-align:left;"> 16:00 </td>
   <td style="text-align:left;"> 20:00 </td>
   <td style="text-align:right;background-color: red !important;"> 31 </td>
   <td style="text-align:right;"> 30 </td>
   <td style="text-align:right;"> 39 </td>
   <td style="text-align:right;"> 21 </td>
   <td style="text-align:right;background-color: green !important;"> 14 </td>
   <td style="text-align:right;background-color: green !important;"> 31 </td>
   <td style="text-align:right;background-color: green !important;"> 0.6186430 </td>
  </tr>
  <tr>
   <td style="text-align:right;"> 6 </td>
   <td style="text-align:left;"> 20:00 </td>
   <td style="text-align:left;"> 24:00 </td>
   <td style="text-align:right;background-color: red !important;"> 14 </td>
   <td style="text-align:right;"> 13 </td>
   <td style="text-align:right;"> 19 </td>
   <td style="text-align:right;"> 7 </td>
   <td style="text-align:right;background-color: green !important;"> 0 </td>
   <td style="text-align:right;background-color: green !important;"> 14 </td>
   <td style="text-align:right;background-color: green !important;"> 0.6751315 </td>
  </tr>
</tbody>
</table>

