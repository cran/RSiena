allEffects <- read.csv('allEffects.csv', as.is=TRUE)
# add default columns needed internally
allEffects$setting <- rep('', nrow(allEffects))
