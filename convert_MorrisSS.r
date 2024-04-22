# After getting results from runMorris_SS.r use
# this script to save results into a .csv format
# for making figures 

# file where Morris SS results are saved
Rdat_fname = "./Results/2024-04-22_MorrisAnalysis_SS.RData"
load(Rdat_fname)

date_to_save <- Sys.Date()
notes <- readline(prompt = "notes for filename: ")

# amt_gut Morris values
var <- 'amt_gut'
save_fname = paste0("./Results/",
                        date_to_save,
                        "_MorrisSS_EE",
                        "_var-", var,
                        "_notes-", notes,
                        ".csv")
# x <- x_amtgut
# temp <- c(x$mu, x$mu.star, x$sigma)
# rownames <- c(names(x$mu))
# colnames <- c("mu", "mu_star", "sigma")
# morvals <- array(temp, dim = c(23, 3), dimnames = list(rownames, colnames))
write.csv(x_amtgut$ee, file = save_fname)

# plas conc
var <- 'plas_conc'
save_fname = paste0("./Results/",
                        date_to_save,
                        "_MorrisSS_EE",
                        "_var-", var,
                        "_notes-", notes,
                        ".csv")
# x <- x_plasconc
# temp <- c(x$mu, x$mu.star, x$sigma)
# morvals <- array(temp, dim = c(23, 3), dimnames = list(rownames, colnames))
write.csv(x_plasconc$ee, file = save_fname)

# inter conc
var <- 'inter_conc'
save_fname = paste0("./Results/",
                        date_to_save,
                        "_MorrisSS_EE",
                        "_var-", var,
                        "_notes-", notes,
                        ".csv")
# x <- x_interconc
# temp <- c(x$mu, x$mu.star, x$sigma)
# morvals <- array(temp, dim = c(23, 3), dimnames = list(rownames, colnames))
write.csv(x_interconc$ee, file = save_fname)

# muscle conc
var <- 'mus_conc'
save_fname = paste0("./Results/",
                        date_to_save,
                        "_MorrisSS_EE",
                        "_var-", var,
                        "_notes-", notes,
                        ".csv")
# x <- x_muscle
# temp <- c(x$mu, x$mu.star, x$sigma)
# morvals <- array(temp, dim = c(23, 3), dimnames = list(rownames, colnames))
write.csv(x_muscle$ee, file = save_fname)