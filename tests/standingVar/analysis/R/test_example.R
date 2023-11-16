A_cpdata <- A.mat(GT_cpdata)
ans.m <- mmer(cbind(Yield,color)~1,
               random=~ vsr(id, Gu=A_cpdata, Gtc=unsm(2)),
               rcov=~ vsr(units, Gtc=unsm(2)), 
               data=DT_cpdata)

ans.m$sigma$`u:id`
