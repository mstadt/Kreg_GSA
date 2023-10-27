# Run each of the single meal experiments

# Load relevant packages
library(deSolve)
library(rootSolve)

# Get relevant functions
source("set_params.r")
source("init_conds.r")
source("varnames.r")
source("mealmod_KClOnly.r")
source("mealmod_MealOnly.r")
source("mealmod_MealKCl.r")

# Set parameter values
params <- set_params()
varnames <- get_varnames()
t0 = 0
tf = (6 + 8) * 60
tvals = seq(t0,tf,1)
#--------------
# Meal + KCl
#--------------
# Solve SS for initial conditions
IC <- init_conds() # baseline ICs


init_guess <- c(amt_gut = IC$amt_gut,
                conc_plas = IC$conc_plas,
                conc_inter = IC$conc_inter,
                conc_muscle = IC$conc_muscle)

ST <- stode(init_guess, time = 0, func = mealmod_MealKCl,
                parms = params)

# Update initial conditions with SS solution
IC <- as.list(ST$y)

# Run the model simulation
out_MealKCl <- as.data.frame(lsoda(unlist(IC[varnames]),
                                times = tvals,
                                func= mealmod_MealKCl,
                                parms = params,
                                rtol = 1e-8,
                                atol = 1e-10
                                )
                            )

#--------------
# Meal Only
#--------------
# Solve SS for initial conditions
IC <- init_conds() # baseline ICs


init_guess <- c(amt_gut = IC$amt_gut,
                conc_plas = IC$conc_plas,
                conc_inter = IC$conc_inter,
                conc_muscle = IC$conc_muscle)

ST <- stode(init_guess, time = 0, func = mealmod_MealOnly,
                parms = params)

# Update initial conditions with SS solution
IC <- as.list(ST$y)

# Run the model simulation
out_MealOnly <- as.data.frame(lsoda(unlist(IC[varnames]),
                                times = tvals,
                                func= mealmod_MealOnly,
                                parms = params,
                                rtol = 1e-8,
                                atol = 1e-10
                                )
                            )

#--------------
# KCl Only
#--------------
# Solve SS for initial conditions
IC <- init_conds() # baseline ICs


init_guess <- c(amt_gut = IC$amt_gut,
                conc_plas = IC$conc_plas,
                conc_inter = IC$conc_inter,
                conc_muscle = IC$conc_muscle)

ST <- stode(init_guess, time = 0, func = mealmod_KClOnly,
                parms = params)

# Update initial conditions with SS solution
IC <- as.list(ST$y)

# Run the model simulation
out_KClOnly <- as.data.frame(lsoda(unlist(IC[varnames]),
                                times = tvals,
                                func= mealmod_KClOnly,
                                parms = params,
                                rtol = 1e-8,
                                atol = 1e-10
                                )
                            )

#----------------
# Save results
#----------------
save_info = as.integer(readline(prompt = 'do you want to save? (0/1) '))
if (save_info) {
    notes = readline(prompt = "notes for filename: ")
    today <- Sys.Date()
    # fname = paste(today, "_MealSim_",
    #                     "_notes-", notes,
    #                     sep = "")
    fMealKCl <- paste("./Results/", today,
                            "_MealSim_",
                            "_type-", "MealKCl",
                            "_notes-", notes,
                            ".csv",
                            sep = "")
    write.csv(out_MealKCl, file = fMealKCl)

    fMealOnly <- paste("./Results/", today,
                            "_MealSim_",
                            "_type-", "MealOnly",
                            "_notes-", notes,
                            ".csv",
                            sep = "")
    write.csv(out_MealOnly, file = fMealOnly)

    fKClOnly <- paste("./Results/", today,
                        "_MealSim_",
                        "_type-", "MealOnly",
                        "_notes-", notes,
                        ".csv",
                        sep = "")
    write.csv(out_KClOnly, file = fKClOnly)

    # Save all the RData
    fname <- paste("./Results/", today,
                        "_MealSims_all",
                        "_notes-", notes,
                        ".RData",
                        sep = "")
    save.image(file = fname) # save details of workspace
}