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



data(ovarian)
ovarian$age40 <- ovarian$age-40
m1 <- coxph(Surv(futime,fustat>0)~age40,data=ovarian)
basehaz(m1,centered=TRUE)
-log(survfit(m1)$surv)
age0 <- mean(ovarian$age40)
plot(survfit(m1,newdata=data.frame(age40=age0)),ylim=c(0.4,1))


library(mets)
m2 <- mets::phreg(Surv(futime,fustat>0)~age40,data=ovarian)
H0 <- predict(m2,surv=TRUE,X=0)
s <- predict(m2,surv=TRUE,X=age0)
points(surv~time,s,col="blue",pch=16)


library(mets)
m2 <- mets::phreg(Surv(futime,fustat>0)~age+factor(rx),data=ovarian)
H0 <- predict(m2,surv=TRUE,X=cbind(56,0))
s <- predict(m2,surv=TRUE,X=cbind(56,0))

s56 <- predict(m2,surv=TRUE,X=cbind(56,1))



plot(m2,col="red",add=TRUE)





d <- data.frame(start=c(1,2,5,2,1,7,3,4,8,8),
               stop=c(2,3,6,7,8,9,9,9,14,17),
               event=c(1,1,1,1,1,1,1,0,0,0),
               x=c(1,0,0,1,0,1,1,1,0,0),
               z=c(1.0,0,2.0,0,3.0,0,4.0,0,5.0,0))
d$S <- with(d,Surv(start,stop,event))
(m1 <- coxph(S~x,data=d,robust=TRUE,ties="breslow"))
(m2 <- coxph(S~x+z,data=d,robust=TRUE,ties="breslow"))
(m2 <- coxph(S~x+z,data=d,robust=TRUE))

survfit(m,newdata=data.frame(x=c(0,1),stop=5))$surv[3,]
