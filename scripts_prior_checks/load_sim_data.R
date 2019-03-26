# Load 2010 tree census and locations for seed traps 1-200

library(readr)
library(dplyr)
library(stringr)

# Dimensions of plot
xmin <- 0
xmax <- 1000
ymin <- 0
ymax <- 500    

# Edge buffer size
buffer <- 20

load("data/bci.full7.rdata")
traps <- read_csv("data/seed_traps.csv") %>%
    filter(TRAP %in% 1:200) %>%
    filter(X >= xmin + buffer, X <= xmax - buffer, 
           Y >= ymin + buffer, Y <= ymax - buffer)
colnames(traps) <- str_to_lower(colnames(traps))

sp_tab <- read_csv("data/sp_tab.csv")

trees <- mutate(bci.full7, sp_code6 = str_to_upper(sp)) %>%
    inner_join(select(sp_tab, sp_code6, rdbh_mm)) %>%
    filter(status == "A", dbh >= rdbh_mm * 2/3) %>%
    select(id = treeID, sp_code6, x = gx, y = gy, dbh_mm = dbh) %>%
    mutate(rba_cm2 = pi * (dbh_mm / 20)^2)

rm(bci.full7)

