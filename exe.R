traindata=ecoli

cb <- function(df, sep="\t", dec=",", max.size=(200*1000)){
  # Copy a data.frame to clipboard
  write.table(df, paste0("clipboard-", formatC(max.size, format="f", digits=0)), sep=sep, row.names=FALSE, dec=dec)
}

#matching.method= Nearest, Cutpoint, Interval
par(mfrow=c(1,5))
origin=algorithm(data=traindata,target="Class",method="origin")
smote=algorithm(data=traindata,target="Class",method="smote")
PS.Cut=algorithm(data=traindata,target="Class",method="PS",matching.method="Cutpoint",point=0.2,interval=qnorm(0.90,0,1))
SP.Cut=algorithm(data=traindata,target="Class",method="SP",matching.method="Cutpoint",point=0.2,interval=qnorm(0.90,0,1))
PSP.Cut=algorithm(data=traindata,target="Class",method="PSP",matching.method="Cutpoint",point=0.2,interval=qnorm(0.90,0,1))

measure.origin=Measure(origin$data,traindata,target="Class")
measure.smote=Measure(smote$data,traindata,target="Class")
measure.PS.c=Measure(PS.Cut$data,traindata,target="Class")
measure.SP.c=Measure(SP.Cut$data,traindata,target="Class")
measure.PSP.c=Measure(PSP.Cut$data,traindata,target="Class")

par(mfrow=c(1,5))
origin=algorithm(data=traindata,target="Class",method="origin")
smote=algorithm(data=traindata,target="Class",method="smote")
PS.Interval=algorithm(data=traindata,target="Class",method="PS",matching.method="Interval",point=0.3,interval=qnorm(0.95,0,1))
SP.Interval=algorithm(data=traindata,target="Class",method="SP",matching.method="Interval",point=0.3,interval=qnorm(0.95,0,1))
PSP.Interval=algorithm(data=traindata,target="Class",method="PSP",matching.method="Interval",point=0.3,interval=qnorm(0.95,0,1))


measure.PS.i=Measure(PS.Interval$data,traindata,target="Class")
measure.SP.i=Measure(SP.Interval$data,traindata,target="Class")
measure.PSP.i=Measure(PSP.Interval$data,traindata,target="Class")


par(mfrow=c(1,5))
origin=algorithm(data=traindata,target="Class",method="origin")
smote=algorithm(data=traindata,target="Class",method="smote")
PS.n=algorithm(data=traindata,target="Class",method="PS",matching.method="Nearest",point=0.3,interval=qnorm(0.95,0,1))
SP.n=algorithm(data=traindata,target="Class",method="SP",matching.method="Nearest",point=0.3,interval=qnorm(0.95,0,1))
PSP.n=algorithm(data=traindata,target="Class",method="PSP",matching.method="Nearest",point=0.3,interval=qnorm(0.95,0,1))


measure.PS.n=Measure(PS.n$data,traindata,target="Class")
measure.SP.n=Measure(SP.n$data,traindata,target="Class")
measure.PSP.n=Measure(PSP.n$data,traindata,target="Class")


##output##
mo=measure.origin$records
ms=measure.smote$records




mps.i=measure.PS.i$records
msp.i=measure.SP.i$records
mpsp.i=measure.PSP.i$records
mps.c=measure.PS.c$records
msp.c=measure.SP.c$records
mpsp.c=measure.PSP.c$records
mps.n=measure.PS.n$records
msp.n=measure.SP.n$records
mpsp.n=measure.PSP.n$records





measure.origin$table
measure.smote$table
measure.PS.i$table
measure.SP.i$table
measure.PSP.i$table




cb(mo)
cb(ms)
cb(mps.c)
cb(msp.c)
cb(mpsp.c)
cb(mps.n)
cb(msp.n)
cb(mpsp.n)
cb(mps.i)
cb(msp.i)
cb(mpsp.i)
table(ecoli$Class)




