# Theoretical Partner Choice Model
# Younginger & Friesen, FEMS Microbiology Letters, 2019

# Basic model of allele dynamics & LD for a host that interacts w multiple symbionts
# Host has 1 locus: choosy or not
# Symbiont has 2 loci: signal locus that host locus interacts with, and benefit locus
# Discrete time system tracking the allele dynamics and LD in the symbiont
# Assume externally regulated populations that are of fixed size

# Include recombination? No: "These loci have been shown to be transferred as an island..." Zhang et al. 2014;
# Kaneko et al. 2000; Suominen et al. 2001; Sullivan et al. 2002; Uchiumi et al. 2004; Parker 2012; 
# Haskett et al. 2016
# Also ignoring finite populations, spatial structure

# Four types of symbiont:
# p1A signal, high benefit
# p1a signal, low benefit
# p0A no signal, high benefit
# p0a no signal, low benefit
# LD = p1A*p0a - p0A*p1a !! depends on allele freqs
# b = fraction of benefit the low benefit strain provides (between 0 and 1)
# s = fraction of benefit the low benefit strain receives due to sanctions (between 0 and 1)
# In principle, plants that do better have more resources for symbionts due to comparative advantage of trade
# Add an exponent for this synergy, zero sets it to one.
# so plant resources for symbionts R = (p_A + b*p_a)^t
# Start with R fixed per plant
# A strains get relatively more when there are a strains bc of sanctions: R*(p_A/(p_A+s*p_a))
# a strains get relatively less: R*(s*p_a/(p_A+s*p_a))

# Signaling strains interact w both C & c plants
# RC=1; Rc=1
# w1A = pC*RC*p1A/(p1A+s*p1a) + pc*Rc*p1A/(p1A+p0A+s*p1a+s*p0a)
# w1a = pC*RC*p1a/(p1A+s*p1a) + pc*Rc*p1a/(p1A+p0A+s*p1a+s*p0a)
# Non-signaling strains only interact w non-choosy c plants
# w0A = pc*Rc*p0A/(p1A+p0A+s*p1a+s*p0a)
# w0a = pc*Rc*p0a/(p1A+p0A+s*p1a+s*p0a)

# Host allele dynamics (pC = freq of choosy allele)
# pC(t+1) = pC*wC / (pC*wC + pc*wc)
# pc(t+1) = pc*wc / (pC*wC + pc*wc)
# Choosy host only interacts w symbiont that has signal, but interacts to the same level
# wC = (p1A + p1a*b)/(p1A + p1a)
# wc = p1A + p0A + p1a*b + p0a*b

make.long<-function(out1){
  my.t<-length(out1[,1])
  time<-rep(1:my.t,times=7)
  var<-rep(names(out1),each=my.t)
  data<-c(out1[,1],out1[,2],out1[,3],out1[,4],out1[,5],out1[,6],out1[,7])
  as.data.frame(cbind(time,var,data))->out2
  out2
}

run.basic.model<-function(my.init,tt=20,b=0.5,s=0.5){
  #tt=20
  dat<-as.data.frame(matrix(ncol=7,nrow=tt))
  names(dat)<-c('pC','pc','p1A','p1a','p0A','p0a','LD')
  # Initialize data
  #dat[1,]<-c(0,1,0.25,0.25,0.25,0.25,NA,NA)
  dat[1,1:6]<-my.init
  # = abs(dat[1,]$p1A*dat[1,]$p0a - dat[1,]$p0A*dat[1,]$p1a) # calc initial LD
  dat$LD[1] = abs(dat[1,]$p1A - (dat[1,]$p1A+dat[1,]$p0A)*(dat[1,]$p1A+dat[1,]$p1a))
  print(paste(c('pC','pc','p1A','p1a','p0A','p0a','LD'),dat[1,]))
  print(b); print(s)
  
  #b=0.5 # relative benefit host gets from bad strain
  #s=0.5 # relative benefit bad strain gets due to sanctions
  RC=1; Rc=1 # host resources don't depend on their symbionts
  
  #my.t=1
  for(my.t in 1:tt){
    pC<-dat$pC[my.t]; pc<-dat$pc[my.t]
    p1A<-dat$p1A[my.t]; p1a<-dat$p1a[my.t]; p0A<-dat$p0A[my.t]; p0a<-dat$p0a[my.t] 
    dat[my.t,]
    #print(p1A + p0A + p1a + p0a)
    #print(pC+pc)
    # host fitness and allele dynamics
    wC = (p1A + p1a*b)/(p1A + p1a)
    wc = p1A + p0A + p1a*b + p0a*b
    #print(wC);print(wc)
    pC_next = pC*wC/(pC*wC+pc*wc)
    pc_next = pc*wc/(pC*wC+pc*wc)
    #print(c(pC_next,pc_next,pC_next+pc_next))
    # symbiont fitness and allele dynamics
    w1A = pC*RC/(p1A+s*p1a) + pc*Rc/(p1A+p0A+s*p1a+s*p0a)
    w1a = pC*RC*s/(p1A+s*p1a) + pc*Rc*s/(p1A+p0A+s*p1a+s*p0a)
    w0A = pc*Rc/(p1A+p0A+s*p1a+s*p0a)
    w0a = pc*Rc*s/(p1A+p0A+s*p1a+s*p0a)
    #print(w1A);print(w1a);print(w0A);print(w0a); # when no choosy hosts, no selection on signal
    p1A_next = p1A*w1A/(p1A*w1A+p1a*w1a+p0A*w0A+p0a*w0a)
    p1a_next = p1a*w1a/(p1A*w1A+p1a*w1a+p0A*w0A+p0a*w0a)
    p0A_next = p0A*w0A/(p1A*w1A+p1a*w1a+p0A*w0A+p0a*w0a)
    p0a_next = p0a*w0a/(p1A*w1A+p1a*w1a+p0A*w0A+p0a*w0a)
    #print(c(p1A_next,p1a_next,p0A_next,p0a_next,p1A_next+p1a_next+p0A_next+p0a_next)) # check
    # put into data vector
    dat[my.t+1,1:6]<-c(pC_next,pc_next,p1A_next,p1a_next,p0A_next,p0a_next)
    dat$LD[my.t+1] = abs(p1A_next - (p1A_next+p0A_next)*(p1A_next+p1a_next))
  }
  ggplot(make.long(dat),aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(data)))) +
    geom_line(aes(color = var)) 
  dat
}


library(ggplot2)
library(reshape2)

# c('pC','pc','p1A','p1a','p0A','p0a','LD')
# Case 1: some C due to drift, medium A, new signal
# no choosy hosts, nothing happens
run.basic.model(c(0,1,0.001,0.0,0.499,0.5),tt=20,b=1,s=1)->my.out
ggplot(make.long(my.out),aes(x=as.numeric(as.character(time)),y=as.numeric(as.character(data)))) +
  geom_line(aes(color = var)) 

## Fig. 2A ##
# just 1% choosy hosts, signal increases and so does choosy allele, generates LD
# LD is transient
# no sanctions needed
run.basic.model(c(0.01,0.99,0.001,0.0,0.470,0.5),tt=200,b=.8,s=1)->my.out

my.out$Generations <- c(1:201)

my.out2 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("LD"),
  variable.name = "Groups", value.name = "Frequency")

my.out3 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("LD","pC","pc","p1A","p1a","p0A","p0a"),
  variable.name = "Groups", value.name = "Frequency")

p1 <- ggplot(my.out3, aes(y = Frequency, x = Generations)) + geom_line(aes(color = Groups), size = 1) + 
theme_classic()  

p1 + scale_colour_manual(values = c("LD" = "black","pC" = "#F8766D","pc" = "#B79F00",
  "p1A" = "#00BA38","p1a" = "#00BFC4","p0A" = "#619CFF","p0a" = "#F564E3"), 
  labels = c("LD","pC: choosy host", "pc: non-choosy host", "p1A: w/ signal, good strain", 
  "p1a: w/ signal, poor strain", "p0A: no signal, good strain", "p0a: no signal, poor strain")) + 
  geom_line(data = my.out2, aes(y = Frequency, x = Generations), size = 1) + xlim(0,50) + 
  theme(legend.position = "none", legend.title = element_blank(), axis.text = element_text(size = 14), 
  legend.text = element_text(size = 14), axis.title = element_text(size = 14))

## Fig. 2B ##
# same as Fig. 2A, only 60% benefit from poor strain
run.basic.model(c(0.01,0.99,0.001,0.0,0.470,0.5),tt=200,b=.6,s=1)->my.out

my.out$Generations <- c(1:201)

my.out2 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("LD"),
  variable.name = "Groups", value.name = "Frequency")

my.out3 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("pC","pc","p1A","p1a","p0A","p0a"),
  variable.name = "Groups", value.name = "Frequency")

p2 <- ggplot(my.out3, aes(y = Frequency, x = Generations)) + geom_line(aes(color = Groups), size = 1) + 
  theme_classic()  

p2 + scale_colour_manual(values = c("LD" = "black","pC" = "#F8766D","pc" = "#B79F00",
  "p1A" = "#00BA38","p1a" = "#00BFC4","p0A" = "#619CFF","p0a" = "#F564E3"), 
  labels = c("LD","pC: choosy host", "pc: non-choosy host", "p1A: w/ signal, good strain", 
  "p1a: w/ signal, poor strain", "p0A: no signal, good strain", "p0a: no signal, poor strain")) + 
  geom_line(data = my.out2, aes(y = Frequency, x = Generations), size = 1) + xlim(0,50) + 
  theme(legend.position = "none", legend.title = element_blank(), axis.text = element_text(size = 14), 
  legend.text = element_text(size = 14), axis.title = element_text(size = 14))

## Fig. 2C ##
# just 1% choosy hosts, if the signal mutation is in the bad strain:
# bad strain increases, choosy allele decreases, generates LD
run.basic.model(c(0.01,0.99,0.00,0.001,0.485,0.499),tt=200,b=.8,s=1)->my.out

my.out$Generations <- c(1:201)

my.out2 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("LD"),
  variable.name = "Groups", value.name = "Frequency")

my.out3 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("pC","pc","p1A","p1a","p0A","p0a"),
  variable.name = "Groups", value.name = "Frequency")

p3 <- ggplot(my.out3, aes(y = Frequency, x = Generations)) + geom_line(aes(color = Groups), size = 1) + 
  theme_classic()  

p3 + scale_colour_manual(values = c("LD" = "black","pC" = "#F8766D","pc" = "#B79F00",
  "p1A" = "#00BA38","p1a" = "#00BFC4","p0A" = "#619CFF","p0a" = "#F564E3"), 
  labels = c("LD","pC: choosy host", "pc: non-choosy host", "p1A: w/ signal, good strain", 
  "p1a: w/ signal, poor strain", "p0A: no signal, good strain", "p0a: no signal, poor strain")) + 
  geom_line(data = my.out2, aes(y = Frequency, x = Generations), size = 1) + xlim(0,50) + 
  theme(legend.position = "none", legend.title = element_blank(), axis.text = element_text(size = 14), 
  legend.text = element_text(size = 14), axis.title = element_text(size = 14))


## Fig. 2D ##
# when signal mutation is in the bad strain,
# sanctions needed to maintain beneficial strain, 
# signal and choosy alleles go extinct
run.basic.model(c(0.01,0.99,0.00,0.001,0.485,0.499),tt=200,b=.8,s=.9)->my.out

my.out$Generations <- c(1:201)

my.out2 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("LD"),
  variable.name = "Groups", value.name = "Frequency")

my.out3 <- melt(my.out, id.vars= ("Generations"), measure.vars = c("pC","pc","p1A","p1a","p0A","p0a"),
  variable.name = "Groups", value.name = "Frequency")

p4 <- ggplot(my.out3, aes(y = Frequency, x = Generations)) + geom_line(aes(color = Groups), size = 1) + 
  theme_classic()  

p4 + scale_colour_manual(values = c("LD" = "black","pC" = "#F8766D","pc" = "#B79F00",
  "p1A" = "#00BA38","p1a" = "#00BFC4","p0A" = "#619CFF","p0a" = "#F564E3"), 
  labels = c("LD","pC: choosy host", "pc: non-choosy host", "p1A: w/ signal, good strain", 
  "p1a: w/ signal, poor strain", "p0A: no signal, good strain", "p0a: no signal, poor strain")) + 
  geom_line(data = my.out2, aes(y = Frequency, x = Generations), size = 1) + xlim(0,50) + 
  theme(legend.position = "none", legend.title = element_blank(), axis.text = element_text(size = 14), 
  legend.text = element_text(size = 14), axis.title = element_text(size = 14))

