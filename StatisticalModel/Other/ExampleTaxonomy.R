library(taxize)
library(worrms)

# If you want to go directly to WoRMS there are a number of ways. 

## Get the classification information
# by ID nunber
aphiaID <- 104464
name <- wm_classification(id = aphiaID)

# getting full classification tree by name
aphiaNAME <- "Calanus finmarchicus"
ID <- wm_classification(name = aphiaNAME)

#Use the underscore function to get a series of data at once.
aphiaID <- c(104464, 346030)
name2 <- wm_classification_(id = aphiaID)
# But you'll notice that it stacks the calssifications 
# so you will need to split by the ID column which 
# is the aphiaID you originally searched for

# This example gives a list back, but also more information. 
name3 <- wm_records_taxamatch(aphiaNAME)
# I find this one harder to use. If you search for 100 taxa, 
# you will end up a list of 100 tibbles, which you can access 
# for example, as: name3[[1]]$genus

# I haven't tried a lot of the wm_* functions so explore and see what you can find in the package.


# If you want more flexibility to use different databases, then you can use taxize, which allows 
# you to define the database. The options you can use are determined by the database you choose. 
# I thinkk it is pretty much a wrapper for worrms.

id2 <- classification(aphiaNAME,db="worms")
name4 <- classification(aphiaID,db="worms") # Note that multiple IDs here also give you a list of classification trees.



