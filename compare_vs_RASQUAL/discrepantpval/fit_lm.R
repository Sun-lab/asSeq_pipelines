#setwd("C:/research/R01/reqtl/2019_05_10")
rasq = read.csv("linear_model_cov.csv")
rasq$rnlt = -log10(rasq$pvalR)#neglog10(rasq)
rasq$tnlt = -log10(rasq$pvalT)#neglog10(trecase)

odNB = rasq$odNB; odNB[odNB<1e-3]=1e-3; odNB[odNB>8]=8; 
odBB = rasq$odBB; odBB[odBB<1e-3]=1e-3; odBB[odBB>8]=8
lodNB = cent(norm(log(odNB))); #lodNB=lodNB-min(lodNB)#-log(min(odNB)))
lodBB = cent(norm(log(odBB))); #lodBB=lodBB-min(lodBB)
lodR = cent(norm(-log(rasq$OD)));



hiR = function(cuts, y){
   kp = abs(y)>cuts
   c(mean(y[kp]>0), sum(kp)) 
}
#show RASQUAL 
GfSNP = rasq$NfSNP
GfSNP[rasq$NfSNP>1]=3
GfSNP[rasq$NfSNP>2]=6
GfSNP[rasq$NfSNP>4]=12
GfSNP[rasq$NfSNP>8]=24
GfSNP[rasq$NfSNP>16]=48
cutoff = 20
rasq$rnlt[rasq$rnlt>cutoff] = cutoff
rasq$tnlt[rasq$tnlt>cutoff] = cutoff
y = rasq$rnlt-rasq$tnlt
z = y>0
kp0 = abs(y)>=0
kp1 = abs(y)>=1
kp5 = abs(y)>=5
kp10 = abs(y)>=10
kp15 = abs(y)>=15

ag0 = aggregate(z[kp0], by=list(GfSNP[kp0]), FUN=mean)
ag1 = aggregate(z[kp1], by=list(GfSNP[kp1]), FUN=mean)
ag5 = aggregate(z[kp5], by=list(GfSNP[kp5]), FUN=mean)
ag10 = aggregate(z[kp10], by=list(GfSNP[kp10]), FUN=mean)
ag15 = aggregate(z[kp15], by=list(GfSNP[kp15]), FUN=mean)

cg0 = aggregate(rep(1,sum(kp0)), by=list(GfSNP[kp0]), FUN=sum)
cg1 = aggregate(rep(1,sum(kp1)), by=list(GfSNP[kp1]), FUN=sum)
cg5 = aggregate(rep(1,sum(kp5)), by=list(GfSNP[kp5]), FUN=sum)
cg10 = aggregate(rep(1,sum(kp10)), by=list(GfSNP[kp10]), FUN=sum)
cg15 = aggregate(rep(1,sum(kp15)), by=list(GfSNP[kp15]), FUN=sum)
ind = 1:15
frR = sapply(ind, hiR, y=y)

numg = cors = cors2 = odnb = rep(0, length(ind))
for(i in ind){
kpi = abs(y)>i#using top 1527 significantly different explains 39% of variation
numg[i] = sum(kpi)
cors[i] = round((cor(lodR[kpi], lodNB[kpi])),2)
cors2[i] = round((cor(lodR[kpi], lodBB[kpi])),2)
odnb[i] = mean(lodNB[kpi])
}

cexes = 2
#png("RASQ_discr_2x2.png", height=8, width=8, res=300, units="in")
#par(mfrow=c(2,2), mar=c(5,5,3,1))

plot(ind, frR[1,], xlab="|log10(T p-val)-log10(R p-val)|", ylab="%p(R)<p(T)", 
main="", bty="n", cex.main=cexes, cex.lab=cexes, cex.axis=cexes, cex=log10(frR[2,]))
legend("topleft", legend=c("#genes", 10^(4:2)), pt.cex=c(1, 4:2), cex=cexes, bty="n", pch=c(NA, 1,1,1))
legend("bottomright", legend=c(2500, 500, 100), pt.cex=log10(c(2500, 500, 100)), cex=cexes, bty="n", pch=c(1,1,1))
axis(1, at=c(12), 12, cex.axis=cexes)

plot(ind, odnb, xlab="|log10(T p-val)-log10(R p-val)|", ylab="normalized log(NBod)", 
main="", bty="n", cex.main=cexes, cex.lab=cexes, cex.axis=cexes, cex=log10(frR[2,]))
legend("bottomleft", legend=c(2500, 500, 100), pt.cex=log10(c(2500, 500, 100)), cex=cexes, bty="n", pch=c(1,1,1))
axis(1, at=c(12), 12, cex.axis=cexes)


x = 1:6
plot(x, ag0[,2], ylim=c(0,1), type="b", xlab="#fSNP", ylab="%p(R)<p(T)", 
main="", bty="n", xaxt="n", cex=log10(cg0[,2]), xlim=c(1,6.1), cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
points(x, ag5[,2], col="blue", type="b", pch=2, cex=log10(cg5[,2]))
points(x, ag10[,2], col="goldenrod", type="b", pch=3, cex=log10(cg10[,2]))
points(x, ag15[,2], col="red", type="b", pch=4, cex=log10(cg15[,2]))
axis(1, at=x, ag0[,1], cex.lab=cexes, cex.axis=cexes)
legend("topleft", legend=c("all",">5", ">10", ">15"), 
text.col=c("black", "blue", "goldenrod", "red"), pch=c(1:4), bty="n", cex=cexes)


plot(ind, cors, cex=log10(numg), pch=1, bty="n", xlab="|log10(T p-val)-log10(R p-val)|", 
ylab="OD correlation", main="", cex.main=cexes, cex.lab=cexes, cex.axis=cexes)
points(ind, cors2, cex=log10(numg), pch=2)
axis(1, at=c(12), 12, cex.axis=cexes)
legend("bottomleft", c("OD(NB)", "OD(BB)"), bty="n", pch=c(1,2), cex=cexes)
legend("topright", legend=c("#genes", 10^(4:2)), pt.cex=c(1, 4:2), cex=cexes, bty="n", pch=c(NA, 1,1,1))
dev.off()




#fit model
cent = function(x){
  med = median(x)
  (x-med)
}
norm = function(x){
  sd = sd(x)
  x/sd
}
odNB = rasq$odNB; odNB[odNB<1e-3]=1e-3; odNB[odNB>8]=8; 
odBB = rasq$odBB; odBB[odBB<1e-3]=1e-3; odBB[odBB>8]=8
odBBf = odBB==1e-3
nodNB = cent(norm(odNB))
nodBB = cent(norm(odBB))
nodNBB = nodNB*nodBB

NfSNP = rasq$NfSNP;
NfSNP[rasq$NfSNP>=60] = 60;
NfSNP = log(NfSNP)
NfSNP = cent(norm(NfSNP))

NrSNP = rasq$NrSNP;
NrSNP = log(NrSNP)
NrSNP = cent(norm(NrSNP))
NrSNP[NrSNP< -5] = -5
NrSNP[NrSNP> +5] =  5
hist(NrSNP)

af = rasq$AF
af = norm(af)

hwec2 = rasq$HWEChisq; hwec2[hwec2<1e-2]=1e-2; hwec2[hwec2>30]=30;hwec2=log(hwec2)
hwec2 = norm(hwec2)

derr = rasq$DeltaErr; derr[derr<1e-10]=1e-10; derr[derr>0.01]=0.01
derr = norm(derr)
phib = abs(rasq$PhiBias-.5); phib[phib>.1]=.1
phib = norm(phib)

asT = rasq$asT
asT = cent(norm(log(asT)))

asR = cent(norm(log(rasq$asR)))

cutoff = 15
rasq$rnlt[rasq$rnlt>cutoff] = cutoff
rasq$tnlt[rasq$tnlt>cutoff] = cutoff
y = rasq$rnlt-rasq$tnlt
kp = abs(y)>=(cutoff-5)

mpp = rasq$mpp
mpp[mpp<1e-5] = 1e-5
mpp = norm(mpp)

#interactions
i12 = asR*NfSNP
i13 = asR*lodBB
i23 = NfSNP*lodBB
i123 = asR*lodBB*NfSNP
ia14 = asT*asR

nms = c("log(OD_BB)", "ASReC_R", "n-fSNP", "ASReC_R:n-SNP", "ASReC_R:log(OD_BB)",
 "n-fSNP:log(OD_BB)", "ASReC_R:n-fSNP:log(OD_BB)", "n-rSNP", "log(OD_NB)",
 "ASReC_T", "ASReC_T:ASReC_R", "AF", "Chi2-HWEC", "Map.Err", "Ref.Bias", "Med(perm-p)") 
lm.2 = lm(y[kp]~lodBB[kp]+asR[kp]+NfSNP[kp]
           +i12[kp]+i13[kp]+i23[kp]+i123[kp]
           +NrSNP[kp]+lodNB[kp]+asT[kp]+ia14[kp]
           +af[kp]+hwec2[kp]+derr[kp]+phib[kp]+mpp[kp])
names(lm.2$coefficients)[-1] = nms
summary(lm.2)



library(car)
an.1 = anova(lm.2)
tss = sum(an.1[,2])
an.3 = Anova(lm.2, type="3")
ano3 = cbind(round(100*an.3[,1]/tss,1), format(an.3[,4], scientific=T, digits=2))[-1,]


lm.m0 = lm(y[kp]~1)
lm.m1 = lm(y[kp]~lodBB[kp])
lm.m2 = lm(y[kp]~asR[kp])
lm.m3 = lm(y[kp]~NfSNP[kp])
lm.m4 = lm(y[kp]~i12[kp])
lm.m5 = lm(y[kp]~i13[kp])
lm.m6 = lm(y[kp]~i23[kp])
lm.m7 = lm(y[kp]~i123[kp])
lm.m8 = lm(y[kp]~NrSNP[kp])
lm.m9 = lm(y[kp]~lodNB[kp])
lm.m10 = lm(y[kp]~asT[kp])
lm.m11 = lm(y[kp]~ia14[kp])

lm.m12 = lm(y[kp]~af[kp])
lm.m13 = lm(y[kp]~hwec2[kp])
lm.m14 = lm(y[kp]~derr[kp])
lm.m15 = lm(y[kp]~phib[kp])
lm.m16 = lm(y[kp]~mpp[kp])

mr2 = 
c(round(summary(lm.m0)$r.squared*100,1), lm.m0$coef[1], summary(lm.m0)$coefficients[1,4],
round(summary(lm.m1)$r.squared*100,1), lm.m1$coef[2], summary(lm.m1)$coefficients[2,4],
round(summary(lm.m2)$r.squared*100,1), lm.m2$coef[2], summary(lm.m2)$coefficients[2,4],
round(summary(lm.m3)$r.squared*100,1), lm.m3$coef[2], summary(lm.m3)$coefficients[2,4],
round(summary(lm.m4)$r.squared*100,1), lm.m4$coef[2], summary(lm.m4)$coefficients[2,4],
round(summary(lm.m5)$r.squared*100,1), lm.m5$coef[2], summary(lm.m5)$coefficients[2,4],
round(summary(lm.m6)$r.squared*100,1), lm.m6$coef[2], summary(lm.m6)$coefficients[2,4],
round(summary(lm.m7)$r.squared*100,1), lm.m7$coef[2], summary(lm.m7)$coefficients[2,4],
round(summary(lm.m8)$r.squared*100,1), lm.m8$coef[2], summary(lm.m8)$coefficients[2,4],
round(summary(lm.m9)$r.squared*100,1), lm.m9$coef[2], summary(lm.m9)$coefficients[2,4],
round(summary(lm.m10)$r.squared*100,1), lm.m10$coef[2], summary(lm.m10)$coefficients[2,4],
round(summary(lm.m11)$r.squared*100,1), lm.m11$coef[2], summary(lm.m11)$coefficients[2,4],
round(summary(lm.m12)$r.squared*100,1), lm.m12$coef[2], summary(lm.m12)$coefficients[2,4],
round(summary(lm.m13)$r.squared*100,1), lm.m13$coef[2], summary(lm.m13)$coefficients[2,4],
round(summary(lm.m14)$r.squared*100,1), lm.m14$coef[2], summary(lm.m14)$coefficients[2,4],
round(summary(lm.m15)$r.squared*100,1), lm.m15$coef[2], summary(lm.m15)$coefficients[2,4],
round(summary(lm.m16)$r.squared*100,1), lm.m16$coef[2], summary(lm.m16)$coefficients[2,4])
mr2 = t(matrix(mr2, nrow=3))

tab = summary(lm.2)$coef
tab = tab[,-3]
tab = cbind(tab, mr2[,2], mr2[,3])
tab[,1] = round(tab[,1],3)
tab[,2] = round(tab[,2],3)
tab[,3] = format(as.numeric(tab[,3]), scientific=T, digits=2)
tab[,4] = round(as.numeric(tab[,4]), 3)
tab[,5] = format(as.numeric(tab[,5]), scientific=T, digits=2)
colnames(tab) = c("Est", "SE", "P-val", "Marg.Est", "Marg.P")
tab

ano = an.1
ano[,3] = round(100*ano[,2]/tss,1)
ano[,5] = format(ano[,5], scientific=T, digits=2)
ano = cbind(ano[,-c(1,4)], ano3)
ano = ano[-17,]
ano = cbind(ano, mr2[-1,1])
rownames(ano) = names(lm.2$coefficients)[-17]
colnames(ano)=c("SumSq", "Typ1per", "T1P-val", "Typ3per", "T3P-val", "Marg.R")
ano

write.csv(ano, "anova.csv", quote=F)
write.csv(tab, "model.csv", quote=F)


