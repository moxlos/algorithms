library(data.table)
library(lubridate)

ABTestDyn <- setRefClass(
  "ABTest",
  fields = list(
    df = "ANY"
  ),
  methods = list(
    initialize = function(df)
    {
      #' df column names:
      #' group : string
      #' timestamp : timestamp
      #' clickedTrue : factor
      df <<- df
    },
    plot_mos = function()
    {
      mosaicplot(table(df[,.(group, clickedTrue)]), main='AB')
    },
    prc_diff = function()
    {
      cvr = df[,.(clickedTrue = mean(as.integer(clickedTrue)-1)), by=.(group)]
      return((cvr$clickedTrue[1] - cvr$clickedTrue[2])/cvr$clickedTrue[1] *100)
    },
    p_value = function(from_dt, to_dt)
    {
      return(prop.test(table(df[timestamp>= from_dt & timestamp <= to_dt, .(group, clickedTrue)]))$p.value)
    }
  )
)


main = function(){
  # Choose parameters:
  pA <- 0.06 # True click through rate for group A
  pB <- 0.06 # True click through rate for group B
  nA <- 500 # Number of cases for group A
  nB <- 500 # Number of cases for group B
  
  
  set.seed(47849)
  df <- data.table(group = rep(c("A", "B"), c(nA, nB)),
                   timestamp = sample(seq(as_datetime('2016-06-02'),
                                          as_datetime('2016-06-09'), by = 1), nA+nB),
                   clickedTrue = as.factor(c(rbinom(n = nA, size = 1, prob = pA),
                                             rbinom(n = nB, size = 1, prob = pB))))
  
  
  
  # Order data by timestamp
  setorder(df, timestamp)
  levels(df$clickedTrue) <- c("0", "1")
  
  print(head(df))
  
  ab <- ABTestDyn$new(df)
  
  ab$plot_mos()
  cat('% difference ',ab$prc_diff())
  cat('significance ', ab$p_value(min(df$timestamp),
                                  '2016-06-07'))
  
  
}


if (getOption('run.main', default=TRUE)) {
  
  main()
}



