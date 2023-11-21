source('code/0_packages.R')
source('code/0_functions.R')

# import datasets ####
train_dataset <- read_excel('data/training_dataset.xlsx')
test_dataset  <- read_excel(     'data/test_dataset.xlsx')

train_dataset$cp <- relevel(factor(train_dataset$cp), ref = 4)
test_dataset$cp  <- relevel(factor(test_dataset$cp),  ref = 4) 

train_dataset$ecg <- relevel(factor(train_dataset$ecg), ref = "0")
test_dataset$ecg  <- relevel(factor(test_dataset$ecg),  ref = "0")

#summary(train_dataset)
#summary(test_dataset)
#
# extract names and indices of new markers ####
ref_vars <- c("sex", "age", "bp", "chol")
new_vars_all <- new.vars(ref_vars)
new_vars_index_all <- new.vars.index(new_vars_all) #42 

# new vars for independent markers
new_vars <- new_vars_all[-(19:36)]
new_vars_index <- new_vars_index_all
new_vars_index[25:42] <- FALSE

# new_vars for correlated markers
new_vars_corr <- new_vars_all[19:36]
new_vars_index_corr <- new_vars_index_all
new_vars_index_corr[1:24] <- FALSE
#
# build ref and new models ####
ref_model <- build.ref.model(ref_vars, train_dataset)
summary(ref_model) 
new_models <- build.new.models(ref_vars, new_vars, train_dataset)
names(new_models) <- new_vars
#
# build ROC curves ####
ROC_ref <- build.ROC.ref(ref_model) 
ROC_new <- build.ROC.new(new_models) 

ROC_ref_test <- build.ROC.ref.test(ref_model,  test_dataset) 
ROC_new_test <- build.ROC.new.test(new_models, test_dataset)
#
# comparisons ####
comparisons_train <- compare.ref.with.new.models(ref_model, new_models, ROC_ref, ROC_new, new_vars)
comparisons_train

# probabilities ####
probs_train <- calculate.probs(ref_model, new_models)
probs_test <- calculate.probs.test(ref_model, new_models, test_dataset)
#
# RB & I ####
results_rbi_train <- calculate.rbi(probs_train)
results_rbi_test  <- calculate.rbi(probs_test)
#
# U-smile plots ####
rbi_for_usmile_train <- prepare.rbi.for.usmile(results_rbi_train, new_vars, new_vars_index, comparisons_train)
usmile_train <- plot.usmile(rbi_for_usmile_train, results_rbi_train, results_rbi_test, "U-smile plots (training dataset)")
grid.newpage()
grid.draw(usmile_train)

rbi_for_usmile_test <- prepare.rbi.for.usmile(results_rbi_test, new_vars, new_vars_index, comparisons_train)
usmile_test <- plot.usmile(rbi_for_usmile_test, results_rbi_train, results_rbi_test, "U-smile plots (test dataset)")
grid.newpage()
grid.draw(usmile_test)
#
# ROC plots ####
plot.roc(ROC_ref, ROC_new, new_vars_index, comparisons_train, title="ROC curves (training dataset)")
plot.roc(ROC_ref_test, ROC_new_test, new_vars_index, comparisons_train, title="ROC curves (test dataset)")
#
# PIW plots ####
plot.piw(probs_train, new_vars_index, comparisons_train, title="PIW plots (training dataset)")
plot.piw(probs_test, new_vars_index, comparisons_train, title="PIW plots (test dataset)")
#
# analysis for correlated markers ####
new_models_corr <- build.new.models(ref_vars, new_vars_corr, train_dataset)
names(new_models_corr) <- new_vars_corr

# build ROC curves 
ROC_new_corr <- build.ROC.new(new_models_corr) 
ROC_new_test_corr <- build.ROC.new.test(new_models_corr, test_dataset)

# comparisons 
comparisons_train_corr <- compare.ref.with.new.models(ref_model, new_models_corr, ROC_ref, ROC_new_corr, new_vars_corr)
comparisons_train_corr

# probabilities
probs_train_corr <- calculate.probs(ref_model, new_models_corr)
probs_test_corr <- calculate.probs.test(ref_model, new_models_corr, test_dataset)

# RB & I
results_rbi_train_corr <- calculate.rbi(probs_train_corr)
results_rbi_test_corr <-  calculate.rbi(probs_test_corr)

# U-smile plots
rbi_for_usmile_train_corr <- prepare.rbi.for.usmile(results_rbi_train_corr, new_vars_corr, new_vars_index_corr, comparisons_train_corr)
usmile_train_corr <- plot.usmile(rbi_for_usmile_train_corr, results_rbi_train_corr, results_rbi_test_corr, "U-smile plots (training dataset)")
#grid.newpage()
#grid.draw(usmile_train_corr)

rbi_for_usmile_test_corr <- prepare.rbi.for.usmile(results_rbi_test_corr, new_vars_corr, new_vars_index_corr, comparisons_train_corr)
usmile_test_corr <- plot.usmile(rbi_for_usmile_test_corr, results_rbi_train_corr, results_rbi_test_corr, "U-smile plots (test dataset)")
#grid.newpage()
#grid.draw(usmile_test_corr)

# combine U-smile plots 
a <- ggarrange(usmile_train_corr, usmile_test_corr, ncol=1, nrow=2)
grid.newpage()
grid.draw(a)

# ROC plots
plot.roc(ROC_ref, ROC_new_corr, new_vars_index_corr, comparisons_train_corr, title="ROC curves (training dataset)")
plot.roc(ROC_ref_test, ROC_new_test_corr, new_vars_index_corr, comparisons_train_corr, title="ROC curves (test dataset)")


# PIW plots 
plot.piw(probs_train_corr, new_vars_index_corr, comparisons_train_corr, title="PIW plots (training dataset)")
plot.piw(probs_test_corr, new_vars_index_corr, comparisons_train_corr, title="PIW plots (test dataset)")



