# After getting results from "runODEMorris.r",
# use this script to save the results into a .csv format
# for making figures

# file where Morris results are saved
Rdat_fname = "./Results/2024-04-22_MorrisAnalysis_MealOnly.RData"
obj_type <- "MealOnly" #"MealOnly"
load(Rdat_fname) # load data into workspace

date_to_save <- Sys.Date()
notes <- readline(prompt = "notes for filename: ")

# amt_gut
var <- "amt_gut"
save_fname = paste("./Results/", 
                    date_to_save,
                    "_MorrisAnalysis",
                    "_type-", obj_type,
                    "_var-", var,
                    "_notes-", notes,
                    ".csv",
                    sep = "")
write.csv(Kmod_res_morris$amt_gut, file = save_fname)

#conc_plas
var <- "conc_plas"
save_fname = paste("./Results/", 
                    date_to_save,
                    "_MorrisAnalysis",
                    "_type-", obj_type,
                    "_var-", var,
                    "_notes-", notes,
                    ".csv",
                    sep = "")
write.csv(Kmod_res_morris$conc_plas, file = save_fname)

# conc_inter
var <- "conc_inter"
save_fname = paste("./Results/", 
                    date_to_save,
                    "_MorrisAnalysis",
                    "_type-", obj_type,
                    "_var-", var,
                    "_notes-", notes,
                    ".csv",
                    sep = "")
write.csv(Kmod_res_morris$conc_inter, file = save_fname)

# conc_muscle
var <- "conc_muscle"
save_fname = paste("./Results/", 
                    date_to_save,
                    "_MorrisAnalysis",
                    "_type-", obj_type,
                    "_var-", var,
                    "_notes-", notes,
                    ".csv",
                    sep = "")
write.csv(Kmod_res_morris$conc_muscle, file = save_fname)