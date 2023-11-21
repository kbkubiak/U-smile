# read my functions ####
source('code/0_packages.R')
#
## create training and test datasets ##
# variable names ####
# original names       new names
# 1. #3  (age)       'age',
# 2. #4  (sex)       'sex',
# 3. #9  (cp)        'cp',
# 4. #10 (trestbps)  'bp',
# 5. #12 (chol)      'chol',
# 6. #16 (fbs)       'glu',
# 7. #19 (restecg)   'ecg',
# 8. #32 (thalach)   'hr',
# 9. #38 (exang)     'exang',
# 10. #40 (oldpeak)  'stde',
# 11. #41 (slope)    'slope',
# 12. #44 (ca)       'ca',
# 13. #51 (thal)     'thal',
# 14. #58 (num)      'diag'            (the predicted attribute)

# import raw data ####
colnames <- c('age',
              'sex',
              'cp',
              'bp',
              'chol',
              'glu',
              'ecg',
              'hr',
              'exang',
              'stde',
              'slope',
              'ca',
              'thal',
              'num',
              'location')

raw_cl <- read.csv("data/raw_data/processed.cleveland.data", header=FALSE, na.strings="?")
raw_hu <- read.csv("data/raw_data/processed.hungarian.data", header=FALSE, na.strings="?")
raw_sw <- read.csv("data/raw_data/processed.switzerland.data", header=FALSE, na.strings="?")
raw_va <- read.csv("data/raw_data/processed.va.data", header=FALSE, na.strings="?")

raw_cl$location <- rep('cl', nrow(raw_cl))
raw_hu$location <- rep('hu', nrow(raw_hu))
raw_sw$location <- rep('sw', nrow(raw_sw))
raw_va$location <- rep('va', nrow(raw_va))

raw_data <- rbind(raw_cl,
                  raw_hu,
                  raw_sw,
                  raw_va)

colnames(raw_data) <- colnames
#summary(raw_data)

disease <- ifelse(raw_data$num == 0, 0, 1)
raw_data <- cbind(disease, raw_data)

sapply(raw_data, function(x) round(sum(is.na(x)) / length(x) * 100, 1)) 
# drop slope, ca, thal - too many NAs: 33.6, 66.4, 52.8 %
# and drop num
raw_data <- select(raw_data, -slope, -ca, -thal, -num)

raw_data <- na.omit(raw_data)
#summary(raw_data) # bp and chol = 0 - must be removed

raw_data <- raw_data[-c(which(raw_data$bp == 0), which(raw_data$chol == 0)),]
summary(raw_data)
#
# generate rnd variables ####
n <- nrow(raw_data)
n0 <- sum(raw_data$disease == 0)
n1 <- sum(raw_data$disease == 1)

# generate independent vars
set.seed(101808)

rnd_normal   <- rnorm(n)
rnd_uniform  <- runif(n, 0, 10)
rnd_exp      <- rexp(n, 1)
rnd_bernoulli <- rbinom(n, 1, 0.8)
rnd_binomial <- rbinom(n, 6, 0.8)
rnd_poisson  <- rpois(n, 1)
strat_rnd_normal   <- ifelse(raw_data$disease == 0, rnorm(n0, 10, 2), rnorm(n1, 12, 2))
strat_rnd_uniform  <- ifelse(raw_data$disease == 0, runif(n0, 0, 6), runif(n1, 2, 8))
strat_rnd_exp      <- ifelse(raw_data$disease == 0, rexp(n0, 0.5), rexp(n1, 1))
strat_rnd_bernoulli <- ifelse(raw_data$disease == 0, rbinom(n0, 1, 0.5), rbinom(n1, 1, 0.2))
strat_rnd_binomial <- ifelse(raw_data$disease == 0, rbinom(n0, 7, 0.6), rbinom(n1, 7, 0.5))
strat_rnd_poisson  <- ifelse(raw_data$disease == 0, rpois(n0, 1), rpois(n1, 1.6))

raw_data_rnd <- cbind(raw_data,
                      rnd_normal,
                      rnd_uniform,
                      rnd_exp,
                      rnd_bernoulli,
                      rnd_binomial,
                      rnd_poisson,
                      strat_rnd_normal,
                      strat_rnd_uniform,
                      strat_rnd_exp,
                      strat_rnd_bernoulli,
                      strat_rnd_binomial,
                      strat_rnd_poisson)

set.seed(5717383)
# correlated rnd vars
raw_data_rnd$rnd_normal0_1 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.1)
raw_data_rnd$rnd_normal0_2 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.2)
raw_data_rnd$rnd_normal0_3 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.3)
raw_data_rnd$rnd_normal0_4 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.4)
raw_data_rnd$rnd_normal0_5 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.5)
raw_data_rnd$rnd_normal0_6 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.6)
raw_data_rnd$rnd_normal0_7 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.7)
raw_data_rnd$rnd_normal0_8 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.8)
raw_data_rnd$rnd_normal0_9 <- rnorm_pre(raw_data_rnd$age, mu = 0, sd = 1, r = 0.9)

# correlated vars for nonevents 
raw_data_rnd %>% filter(disease == 0) %>% 
  mutate(strat_rnd_normal0_1=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.1),
         strat_rnd_normal0_2=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.2),
         strat_rnd_normal0_3=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.3),
         strat_rnd_normal0_4=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.4),
         strat_rnd_normal0_5=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.5),
         strat_rnd_normal0_6=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.6),
         strat_rnd_normal0_7=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.7),
         strat_rnd_normal0_8=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.8),
         strat_rnd_normal0_9=rnorm_pre(age, mu = 10, sd = 2.5, r = 0.9)) -> raw_data_rnd_nonev


# correlated vars for events
raw_data_rnd %>% filter(disease == 1) %>% 
  mutate(strat_rnd_normal0_1=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.1),
         strat_rnd_normal0_2=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.2),
         strat_rnd_normal0_3=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.3),
         strat_rnd_normal0_4=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.4),
         strat_rnd_normal0_5=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.5),
         strat_rnd_normal0_6=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.6),
         strat_rnd_normal0_7=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.7),
         strat_rnd_normal0_8=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.8),
         strat_rnd_normal0_9=rnorm_pre(age, mu = 11, sd = 2.5, r = 0.9)) -> raw_data_rnd_event

raw_data_rnd <- rbind(raw_data_rnd_nonev,
                      raw_data_rnd_event)
summary(raw_data_rnd)

write_xlsx(raw_data_rnd, 'data/raw_data_rnd.xlsx')

# divide into training and test datasets ####
table(raw_data$location,
      raw_data$disease)
train_table <- round(table(raw_data$location,
                           raw_data$disease) * 0.5) # number of non-events and event for the training dataset
cl_0 <- train_table[1,1] #  82
cl_1 <- train_table[1,2] #  70
hu_0 <- train_table[2,1] #  82
hu_1 <- train_table[2,2] #  49
va_0 <- train_table[3,1] #  10
va_1 <- train_table[3,2] #  38

sum(train_table[,1]) # 174 non-events
sum(train_table[,2]) # 157 events

set.seed(101909)

ind_cl_0 <- sample(which(raw_data_rnd$location == 'cl' & raw_data_rnd$disease == 0), cl_0)
ind_cl_1 <- sample(which(raw_data_rnd$location == 'cl' & raw_data_rnd$disease == 1), cl_1)
ind_hu_0 <- sample(which(raw_data_rnd$location == 'hu' & raw_data_rnd$disease == 0), hu_0)
ind_hu_1 <- sample(which(raw_data_rnd$location == 'hu' & raw_data_rnd$disease == 1), hu_1)
ind_va_0 <- sample(which(raw_data_rnd$location == 'va' & raw_data_rnd$disease == 0), va_0)
ind_va_1 <- sample(which(raw_data_rnd$location == 'va' & raw_data_rnd$disease == 1), va_1)

raw_data_rnd[c(ind_cl_0,
               ind_cl_1,
               ind_hu_0,
               ind_hu_1,
               ind_va_0,
               ind_va_1), ] -> raw_train_ds

raw_data_rnd[-c(ind_cl_0,
                ind_cl_1,
                ind_hu_0,
                ind_hu_1,
                ind_va_0,
                ind_va_1), ] -> raw_test_ds

summary(raw_train_ds)
summary(raw_test_ds)

write_xlsx(raw_train_ds, 'data/training_dataset.xlsx')
write_xlsx(raw_test_ds,      'data/test_dataset.xlsx')

