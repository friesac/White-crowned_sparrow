# R language
# You should be able to copy and paste these lines with some minor editing where noted.
# Anything after the "#" will not be read by R

library(adegenet) # use the command install.packages("adegenet") if you've never installed the library before. Then follow up with library(adegenet) which will load it into the workspace so you can call adegenet functions.
library(ggplot)
library(PopGenReport)
library(geosphere)

# You'll change to your own directory.
# setwd("C:/Users/fries/Dropbox/Lisa_Ade/project_data.stru")

# So the format Lisa sent is a basic excel doc which I'm sure there is a way to import into adegenet but I've always used structure files to import into adegenet. Just a change in the loci name row. I've attached what we were working with.
# 285 genotypes, 14 markers, 1 genotype labels, 2 pop factor, blank enter, 1 marker name row, y single row genotypes.
geno_genind <- read.structure("C:/Users/fries/Dropbox/Lisa_Ade/project_data.stru", ask=F, n.ind=285,n.loc=14,onerowperind=T, col.pop=2, col.lab=1, row.marknames=1) 

# Now we'll subset out the unknown individuals from the genind object to test for a.optima and xval functions on the set of KNOWN individuals
# This is based on the file with the first 40 samples being "unknown"

unk_genind <- geno_genind[geno_genind$pop=="P1"] # These are your test set - as long as unknown is in the first row of your file ... P1 will be whatever the first rows pop is.
no_unk_genind <- geno_genind[!geno_genind$pop=="P1"] # These are your training set


#######################################################################################################################
#######################################################################################################################
# Using the a-score to determine The trade-off between power of discrimination and over-fitting

wcsp_dapc <- dapc(no_unk_genind, n.da=100, n.pca=300) # - Looks like 106 is the best number of pcs
optim.a.score(wcsp_dapc)


# Using the xvalDapc cross-validation data - Looks like 100 is the best pc number

mat <-as.matrix(na.replace(no_unk_genind,method="mean"))
grp <- pop(no_unk_genind)

xval <- xvalDapc(mat, grp, n.pca.max=300, training.set =0.90, result="groupMean", center = TRUE, scale=FALSE, n.pca = NULL, n.rep=100, xval.plot=TRUE)

### Both of these solutions suggest that we've gone beyond the n/3 benchmark and are at risk of overfitting... however, they are recommending these.
#######################################################################################################################
#######################################################################################################################

# Now we'll run the dapc using the only North/South to define the disc. functions.
dapc_no_unk_genind <- dapc(no_unk_genind, grp=truenames(no_unk_genind@pop))
# You'll be prompted to pick the number of PCs to retain... looks like around 90% of your total dataset variation is defined in ~75 PCs... don't go higher than N/3 but you could play around with other pc numbers... there are alternate ways to pick PCs but probably won't make much difference.
# Keep the 1 discriminant function

scatter(dapc_no_unk_genind, legend=T) # This is the plot we originally looked at that defined North and South on one axis.

# Now we'll take the function that splits North and South and use it to predict assignment of unknowns
predicted_unknown_assignments <- predict.dapc(dapc_no_unk_genind, newdata=unk_genind)

# The predicted_unknown_assignments is a special R object ("list") that has multiple variables embedded in it... access them by using the "$" sign after the predicted_unknown_assignments variable name... e.g. 

predicted_unknown_assignments$posterior 

# or

predicted_unknown_assignments$ind.score

# will give you the posterior probabilities of assignment for all the unknown individuals to each group OR the actual loading score of the unknown indivduals.

# For the leave one out analysis, in order to test the stability of your original discriminant function to distinguish North and South try this...

how_many_to_leave_out <- 10
iterations <- 10
pcas <- 82
predict_prior <- NULL
predict_post <- NULL
predict_name <- NULL

for (i in 1:iterations) {
	extractor <- sample(1:length(no_unk_genind$ind.names), how_many_to_leave_out)
	leave_out_genind <- no_unk_genind[-extractor]
	predict_left_out <- no_unk_genind[extractor]
	dapc_leave_out_genind <- dapc(leave_out_genind, grp=truenames(leave_out_genind@pop), n.pca=pcas, n.da=1)
	predict_by_dapc <- predict.dapc(dapc_leave_out_genind, newdata=predict_left_out)
	predict_prior <- c(predict_prior, as.character(truenames(predict_left_out)$pop))
	predict_post <- c(predict_post, as.character(predict_by_dapc$assign))
	predict_name <- c(predict_name, predict_left_out$ind.names)
}
predict_leave_out_matrix <- cbind(predict_name, predict_prior, predict_post, predict_prior==predict_post) 


# Another way that only goes through leaving one out in order. 

pcas <- 75
predict_prior <- NULL
predict_post <- NULL
predict_name <- NULL

for (i in 1:length(no_unk_genind$ind.names)) {
	leave_1_out_genind <- no_unk_genind[-i]
	dapc_leave_1_out_genind <- dapc(leave_1_out_genind, grp=truenames(leave_1_out_genind@pop), n.pca=pcas, n.da=1)
	predict_leave_1_out <- predict.dapc(dapc_leave_1_out_genind, newdata=no_unk_genind[i])
	predict_prior <- c(predict_prior, as.character(truenames(no_unk_genind)$pop[i]))
	predict_post <- c(predict_post, as.character(predict_leave_1_out$assign))
	predict_name <- c(predict_name, no_unk_genind[i]$ind.names)	
}
predict_leave_1_out_matrix <- cbind(predict_name,predict_prior,predict_post,predict_prior==predict_post)



### IBD Tests
### For populations
# geno_genpop<-genind2genpop(no_unk_genind)
# Dgen <- dist.genpop(geno_genpop,method=2)
# wcsp_distance <- read.table("C:/Users/fries/Dropbox/Lisa_Ade/mapsites.dat", header=TRUE, sep="\t")
# wcsp_coords <- cbind(wcsp_distance$long, wcsp_distance$lat)
# Dgeo <- dist(wcsp_coords)
# ibd <- mantel.randtest(Dgen,Dgeo)

### For individuals
Dgen <- gd.kosman(no_unk_genind)
Dgeo <- 

mantel.randtest




