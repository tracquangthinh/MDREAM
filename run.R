# Input
# train_data = readRDS("data/train_data")
# train data is a R list including 3 objects: gex (gene expression), drug_response and mut (mutation)
# All information are available from the BeatAMl publication
# Here we show some examples to illustrate the structure of train_data
# > dim(train_data$gex)
# [1] 22521   461
# > train_data$gex[1:5, 1:5]
#                   Patient1    Patient2  Patient3  Patient4  Patient5
# ENSG00000121410 -1.6992571  2.33217444  2.521210 2.2471926 2.2647500
# ENSG00000268895  0.3268287  1.33082084  3.255150 2.7548006 2.1711128
# ENSG00000148584  0.8062588  1.77261652  2.348147 2.0032284 2.1909398
# ENSG00000175899 -0.7175547  0.03904027 -0.760792 0.2802349 4.4788898
# ENSG00000245105  0.5907237 -2.23178024  2.784988 1.4912317 0.4602694
#
# > dim(train_data$drug_response)
# [1] 32263     4
# > head(train_data$drug_response, 3)
#             inhibitor   lab_id       ic50      auc
# 17-AAG (Tanespimycin) Patient1 10.0000000 225.9180
# 17-AAG (Tanespimycin) Patient2  2.7228447 164.5612
# 17-AAG (Tanespimycin) Patient3  0.1165822 107.7196
#
# > dim(train_data$mut)
# [1] 302 461
# > train_data$mut[1:3, 1:3]
#       Patient1 Patient2 Patient3
# A1CF         0        0        1
# A2M          1        0        0
# AARS2        0        0        0

library(parallel)
library(e1071)

source("Rsource.R")
source("generate_pas.R")
source("train.R")

# train_data is publicly available from BeatAML paper
# train_data = readRDS("data/train_data")

gene_set = readRDS("data/gene_set")
pas = generate_pas(train_data, gene_set, n_cores = 1)

models = train(train_data, pas, ic50 = FALSE, 
                  method_1 = "svmRadial", method_2 = "svmRadial", n_cores = 1)
