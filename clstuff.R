cd ~/
chmod 400 yobbo.pem
insts=$(ec2din | grep INSTANCE | awk {'print $4'})

ec2din | grep INSTANCE | awk {'print $4'}
echo $insts | fmt -1 | head -33
s3cmd get s3://disc-output ~/Dropbox/disc\ mod\ backup/Lancet\ MS/Couples\ Model\ Revision\ 121013/Output/main/ -r --skip-existing

ssh -o "StrictHostKeyChecking no" -i  yobbo.pem ubuntu@

92 hom
89 het

# Request 1 spot instance on c1xlarge
ec2rsi ami-8174ffe8 -p 0.70 -k yobbo -n 1 -t c1.xlarge

## sudo killall -u ubuntu

      
## Check out current instances
system("ec2din | grep INSTANCE | awk {'print $4'} > ~/tempec2s.csv")
ints <- read.csv("~/tempec2s.csv", header = F)
ints[,1] <- levels(ints[,1])[ints[,1]]
ints


to.do <- 4
ints[to.do,]

## kill all
cmds <- NULL
for(ii in to.do)
  {
    cmd <- paste("ssh -o \"StrictHostKeyChecking no\" -i  yobbo.pem ubuntu@",
                 ints[ii,],
                 " sudo killall -u ubuntu", sep = "")
    cmds <- c(cmds,cmd)
  }
write(cmds, "~/est.csv")


cmds <- NULL
for(ii in to.do)
  {
    cmd <- paste("ssh -o \"StrictHostKeyChecking no\" -i  yobbo.pem ubuntu@",
                 ints[ii,],
                 " \"cd files; nohup R CMD BATCH '--args group.ind=9 short.test=F all.cores=T num.cores=8 d.nburn=300 d.nthin=3 d.niter=1300 nburn=300 nthin=1 niter=4300 survive=T tell=100 adapt=T term.on.finish=T simul=T K.sim=2000 heter=F' areldisc6.R > foo.out 2> foo.err < /dev/null \"",
                 " &", sep = "")
    cmds <- c(cmds,cmd)
  }
write(cmds, "~/est.csv")



cmds <- NULL
group.ind <- 10
for(ii in to.do)
  {
    cmd <- paste("ssh -o \"StrictHostKeyChecking no\" -i  yobbo.pem ubuntu@",
                 ints[ii,],
                 " \"cd files; nohup R CMD BATCH '--args group.ind=", group.ind,  " short.test=F all.cores=T num.cores=8 d.nburn=300 d.nthin=3 d.niter=1300 nburn=300 nthin=1 niter=8300 survive=T tell=100 adapt=T term.on.finish=T simul=F K.sim=1000 heter=F low.coverage.arv=F partner.arv=F fsd.sens=F' areldisc6.R > foo.out 2> foo.err < /dev/null \"",
                 " &", sep = "")
    cmds <- c(cmds,cmd)
    group.ind <- group.ind + 1
  }
write(cmds, "~/est.csv")


## kill all instances
## for i in `ec2din | grep running | cut -f2`; do ec2kill $i; done

## ec2din | grep INSTANCE | awk {'print $4'}


## "cd files; nohup R CMD BATCH '--args group.ind=9 short.test=F all.cores=T num.cores=8 d.nburn=300 d.nthin=3 d.niter=1000 nburn=400 nthin=1 niter=4400 survive=T tell=100 adapt=T term.on.finish=T simul=T K.sim=2000 heter=T' areldisc6.R > foo.out 2> foo.err < /dev/null" &
