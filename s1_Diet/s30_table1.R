rm(list = ls())
# Set your working directory and data path here

library(data.table)

df <- fread('data/Y_total.csv')
colnames(df)
library(tableone)
myVars <- c('Age', 'Sex', 'Ethnic_group', 'education_level', 'TDI', 'BMI', 'Smoke',
            'IPAQ', 'energy', 'Vitamin', 'fast_time')
catVars <- c('Sex', 'Ethnic_group', 'education_level', 'Smoke', 'IPAQ')


table <- CreateTableOne(vars = myVars,
                        factorVars = catVars,
                        strata = "death",
                        data = df,
                        addOverall = TRUE
);table

table1<- print(table,
               catDigits = 2,contDigits = 2,pDigits = 3,

               showAllLevels=TRUE,
               quote = FALSE,
               noSpaces = TRUE,
               printToggle = TRUE)
write.csv(table1, file = "results/s0_Association/s0_table1.csv")
