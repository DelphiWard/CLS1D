## ---- echo=FALSE---------------------------------------------------------
library(CLS1D)

## ---- out.width = "600px", echo=FALSE------------------------------------
knitr::include_graphics("TC2_landscapes.png")

## ------------------------------------------------------------------------
str(regime1_sw)

## ------------------------------------------------------------------------
x <- ifelse(regime1_sw == "2", 1, 0) #species 2

## ------------------------------------------------------------------------
set.seed(2)
ex <- CLS1D(x, Lhalf=20, d=3, metric="errX", n.samples=20)

## ---- fig.width=6, fig.height=6------------------------------------------
plotCLS(ex, metric="errX", ylim=c(0,1), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(regime2_sw == "2", 1, 0) #species 2 from regime 2
ex <- CLS1D(x, Lhalf=20, d=3, metric="errX", n.samples=20)
plotCLS(ex, metric="errX", ylim=c(0,1), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
x <- floor(runif(6000, 0,5)) #randomly distributed species
ex <- CLS1D(ifelse(x == "2", 1, 0), Lhalf=20, d=3, metric="errX", n.samples=20)
plotCLS(ex, metric="errX", ylim=c(0,1), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(regime2_sw == "2", 1, 0) #species 2 from regime 2
ex <- CLS1D(x, Lhalf=50, d=3, metric="errX", n.samples=20)
plotCLS(ex, metric="errX", ylim=c(0,1.3), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2007 == "CF", 1, 0) #foliose coral from 2007
ex <- CLS1D(x, Lhalf=20, d=3, metric="errX", n.samples=20)
plotCLS(ex, metric="errX", ylim=c(0,1), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2007 == "CF", 1, 0) #foliose coral from 2007
ex <- CLS1D(x, Lhalf=50, d=3, metric="errX", n.samples=20) #max window size = 100cm
plotCLS(ex, metric="errX", ylim=c(0,1.5), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2010 == "CF", 1, 0) #foliose coral from 2010
ex <- CLS1D(x, Lhalf=50, d=3, metric="errX", n.samples=20) #max window size = 100cm
plotCLS(ex, metric="errX", ylim=c(0,1.5), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2010 == "CF", 1, 0) #foliose coral from 2010
ex <- CLS1D(x, Lhalf=50, d=4, metric="errX", n.samples=20) #d = 4
plotCLS(ex, metric="errX", ylim=c(0,1.5), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2010 == "CF", 1, 0) #foliose coral from 2010
ex <- CLS1D(x, Lhalf=50, d=2, metric="errX", n.samples=20) #d = 2
plotCLS(ex, metric="errX", ylim=c(0,1.5), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2010 == "CF", 1, 0) #foliose coral from 2010
ex <- CLS1D(x, Lhalf=50, d=3, metric="errX", n.samples=20, k=10) #k = 10
plotCLS(ex, metric="errX", ylim=c(0,1.5), xticks = 5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- ifelse(bb2010 == "CF", 1, 0) #foliose coral from 2010
ex <- CLS1D(x, Lhalf=50, d=3, metric="errX", n.samples=20, k=2) # k = 2
plotCLS(ex, metric="errX", ylim=c(0,1.5), xticks = 5)

## ------------------------------------------------------------------------
head(Sp2R1_sts)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- CLS1Dsts(Sp2R1_sts, Lhalf=20, metric="errX", n.samples=20)
plotCLS(x, metric="errX", ylim=c(0,1), xticks=5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- CLS1Dsts(Sp2R1_sts, Lhalf=40, metric="errX", n.samples=20)
plotCLS(x, metric="errX", ylim=c(0,1), xticks=5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x<- Sp2R1_sts[,1:3]
x<- CLS1Dsts(x, Lhalf=20, metric="errX", n.samples=20)
plotCLS(x, metric="errX", ylim=c(0,1), xticks=5)

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- CLS1Dsts(Sp2R1_sts, Lhalf=20, metric="PRSq", n.samples=20)
plotCLS(x, metric="PRSq", ylim=c(0,1))

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- CLS1Dsts(Sp2R1_sts, Lhalf=40, metric="PRSq", n.samples=20)
plotCLS(x, metric="PRSq", ylim=c(0,1))

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- CLS1Dsts(Sp2R1_sts, Lhalf=100, metric="PRSq", n.samples=20)
plotCLS(x, metric="PRSq", ylim=c(0,1))

## ---- fig.width=6, fig.height=6------------------------------------------
set.seed(2)
x <- CLS1Dsts(Sp2R2_sts, Lhalf=20, metric="errX", n.samples=20)
plotCLS(x, metric="errX", ylim=c(0,1), xticks=5)

