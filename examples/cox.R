library(mets)
library(survival)

test2 <- list(start=c(1,2,5,2,1,7,3,4,8,8), 
              stop=c(2,3,6,7,8,9,9,9,14,17), 
              event=c(1,1,1,1,1,1,1,0,0,0), 
              x=c(1,0,0,1,0,1,1,1,0,0)) 
test2$stop <- test2$stop+runif(length(test2$stop))*1e-4

c1 <- coxph(Surv(stop, event) ~ x, test2)
c2 <- phreg(Surv(stop, event) ~ x, test2)
c3 <- cox.aalen(Surv(stop, event) ~ prop(x), test2)

c1 <- coxph(Surv(start,stop, event) ~ x, test2)
c2 <- phreg(Surv(start,stop, event) ~ x, test2)
c3 <- cox.aalen(Surv(start,stop, event) ~ prop(x), test2)




predict(c2,surv=TRUE,newdata=data.frame(x=1))
survfit(c1)$surv
diag(survfit(c1,newdata=data.frame(x=test2$x))$surv)
predict(c2,surv=TRUE)

basehaz(c1,centered=FALSE)
predict(c2,surv=FALSE)
basehaz(c1)
