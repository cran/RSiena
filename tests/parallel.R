library(RSiena)
print(packageDescription("RSiena",fields="Repository/R-Forge/Revision"))

##test3
mynet1 <- sienaNet(array(c(tmp3, tmp4),dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(findiff=TRUE, fn = simstats0c, projname='test3',
                       cond=FALSE, nsub=2, n3=100)
print('test3')
system.time(ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE))#,dll='../siena/src/RSiena.dll')
##test4
mymodel$projname <- 'test4'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
print('test4')
system.time(ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE))#,dll='../siena/src/RSiena.dll')
##test7
mynet1 <- sienaNet(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(fn = simstats0c, projname='test7', nsub=2, n3=100,
                       cond=FALSE)
print('test7')
system.time(ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE))#,dll='../siena/src/RSiena.dll')
##test8
mymodel$projname <- 'test8'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
print('test8')
system.time(ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE))#,dll='../siena/src/RSiena.dll')

mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a,type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
myeff$initialValue[94] <- 0.34699930338 ## siena3 starting values differ
##test10
print('test10')
mymodel$projname <- 'test10'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
system.time(ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, verbose=TRUE))
##test11
print('test11')
system.time(data501 <- sienaDataCreateFromSession("s50.csv", modelName="s50"))
system.time(data501e <- sienaDataCreateFromSession("s50e.csv", modelName="s50e"))
system.time(data501paj <- sienaDataCreateFromSession("s50paj.csv", modelName="s50paj"))

model501e <- model.create( projname="s50e", cond=FALSE, nsub=2, n3=100 )
system.time(ans501e <- siena07(model501e, data=data501e$mydata, effects=data501e$myeff,
                   parallelTesting=TRUE, batch=TRUE, verbose=TRUE))
## compare with outputs in parallelchecked/
