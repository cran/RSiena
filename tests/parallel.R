library(RSiena)
print(packageDescription("RSiena",fields="Repository/R-Forge/Revision"))

##test3
mynet1 <- sienaNet(array(c(tmp3, tmp4),dim=c(32, 32, 2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(findiff=TRUE, fn = simstats0c, projname='test3',
                       cond=FALSE, nsub=2, n3=100)
print('test3')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE, parallelTesting=TRUE, silent=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test4
mymodel$projname <- 'test4'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
print('test4')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test7
mynet1 <- sienaNet(array(c(tmp3,tmp4),dim=c(32,32,2)))
mydata <- sienaDataCreate(mynet1)
myeff<- getEffects(mydata)
mymodel<- model.create(fn = simstats0c, projname='test7', nsub=2, n3=100,
                       cond=FALSE)
print('test7')
ans<- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test8
mymodel$projname <- 'test8'
mymodel$cconditional <- TRUE
mymodel$condvarno <- 1
print('test8')
ans <- siena07(mymodel, data=mydata, effects=myeff,  batch=TRUE,
              parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)#,dll='../siena/src/RSiena.dll')
ans
##test9

mynet1 <- sienaNet(array(c(s501, s502, s503), dim=c(50, 50, 3)))
mynet2 <- sienaNet(s50a,type='behavior')
mydata <- sienaDataCreate(mynet1, mynet2)
myeff <- getEffects(mydata)
myeff$initialValue[98] <- 0.34699930338 ## siena3 starting values differ
##test10
print('test10')
mymodel$projname <- 'test10'
mymodel$cconditional <- TRUE
mymodel$condvarno<- 1
ans <- siena07(mymodel, data=mydata, effects=myeff, batch=TRUE,
               parallelTesting=TRUE, silent=TRUE)
##, verbose=TRUE)
ans
##test11
print('test11')
data501 <- sienaDataCreateFromSession("s50.csv", modelName="s50")
data501e <- sienaDataCreateFromSession("s50e.csv", modelName="s50e")
data501paj <- sienaDataCreateFromSession("s50paj.csv", modelName="s50paj")

model501e <- model.create( projname="s50e", cond=FALSE, nsub=2, n3=100 )
ans501e <- siena07(model501e, data=data501e$mydata, effects=data501e$myeff,
                   parallelTesting=TRUE, batch=TRUE, silent=TRUE)
##, verbose=TRUE)
ans501e
## compare with outputs in parallelchecked/
