# Lapiedra 2018 rec 319
# Phenotypic extractions

# Load libraries
library(here)

# Set wd====
dir<-here()

# Females====
L2018F<-read.csv("Data/Lapiedra2018Fem.csv")

###control
L2018F.C<-subset(L2018F, experimental_treatment=="control")

cor.test(L2018F.C$survival_4months, L2018F.C$leioc_exposed_time)
plot(L2018F.C$survival_4months, L2018F.C$leioc_exposed_time)
abline(lm(L2018F.C$leioc_exposed_time~L2018F.C$survival_4months))
#r=0.01779847

cor.test(L2018F.C$survival_4months, L2018F.C$leioc_head_out_31)
plot(L2018F.C$survival_4months, L2018F.C$leioc_head_out_31)
abline(lm(L2018F.C$leioc_head_out_31~L2018F.C$survival_4months))
#r=-0.2503607 but switch sign to reflect boldness

###predator treatment
L2018F.P<-subset(L2018F, experimental_treatment=="predation")

cor.test(L2018F.P$survival_4months, L2018F.P$leioc_head_out_31)
plot(L2018F.P$survival_4months, L2018F.P$leioc_head_out_31)
abline(lm(L2018F.P$leioc_head_out_31~L2018F.P$survival_4months))
#r=0.08023263 switch sign to reflect boldness

cor.test(L2018F.P$survival_4months, L2018F.P$leioc_exposed_time)
plot(L2018F.P$survival_4months, L2018F.P$leioc_exposed_time)
abline(lm(L2018F.P$leioc_exposed_time~L2018F.P$survival_4months))
#r=-0.3005145

# Males====
L2018M<-read.csv("Data/Lapiedra2018Male.csv")

##control
L2018M.C<-subset(L2018M, experimental_treatment=="control")

cor.test(L2018M.C$survival_4months, L2018M.C$leioc_head_out_31)
plot(L2018M.C$survival_4months, L2018M.C$leioc_head_out_31)
abline(lm(L2018M.C$leioc_head_out_31~L2018M.C$survival_4months))
#r=-0.163642 switch sign to reflect boldness

cor.test(L2018M.C$survival_4months, L2018M.C$leioc_exposed_time)
plot(L2018M.C$survival_4months, L2018M.C$leioc_exposed_time)
abline(lm(L2018M.C$leioc_exposed_time~L2018M.C$survival_4months))
#r=0.0197477

###predation
L2018M.P<-subset(L2018M, experimental_treatment=="predation")

cor.test(L2018M.P$survival_4months, L2018M.P$leioc_head_out_31)
plot(L2018M.P$survival_4months, L2018M.P$leioc_head_out_31)
abline(lm(L2018M.P$leioc_head_out_31~L2018M.P$survival_4months))
#r=-0.2393831 switch sign to reflect boldness

cor.test(L2018M.P$survival_4months, L2018M.P$leioc_exposed_time)
plot(L2018M.P$survival_4months, L2018M.P$leioc_exposed_time)
abline(lm(L2018M.P$leioc_exposed_time~L2018M.P$survival_4months))
#r=-0.2489085