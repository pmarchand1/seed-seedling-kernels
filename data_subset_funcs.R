# Functions to extract the data on adult trees, seeds, recruits and seedlings
#  and format it for input to the models

library(dplyr)
library(stringr)


# Load data and define constants ------------------------------------------

# Read in reformatted data files
trees <- read.csv("data/tree_census_intervals.csv", stringsAsFactors = FALSE)
seedlings <- read.csv("data/seedlings_tidy.csv", stringsAsFactors = FALSE)
seeds <- read.csv("data/seeds_merged_calyr.csv", stringsAsFactors = FALSE)
recruits <- read.csv("data/recruits_merged.csv", stringsAsFactors = FALSE)

# Read species-level data
sp_list <- read.csv("sp_tab.csv", stringsAsFactors = FALSE)

# Dimensions of plot
xmin <- 0
xmax <- 1000
ymin <- 0
ymax <- 500    

# Tree census years
tree_years <- c(1982, 1985, 1990, 1995, 2000, 2005, 2010, 2015) 

# Seedling survey years 
seedling_years <-  c(2002:2004, 2006, 2008:2009, 2011:2013)
seedling_dt <- c(1, diff(seedling_years)) # years between surveys


# Helper functions --------------------------------------------------------

# Interpolation functions for tree basal area and survival at year
#  given values at prev_yr and next_yr
ba_interp <- function(prev_ba, next_ba, prev_yr, next_yr, year) {
    # If non-zero values for both years, interpolate assuming geometric growth
    ifelse(prev_ba == 0 | next_ba == 0, pmax(prev_ba, next_ba),
           prev_ba * (next_ba / prev_ba)^((year - prev_yr) / (next_yr - prev_yr)))
}

surv_interp <- function(next_status, prev_yr, next_yr, year) {
    # prev_status is always alive
    # if next_status is dead, prob. of survival linearly decreases over time
    ifelse(next_status == "A", 1,
           runif(length(next_status)) > (year - prev_yr) / (next_yr - prev_yr))
}


# Data subsetting functions -----------------------------------------------

# Subset "trees" data frame for a given species and year
subset_trees <- function(spc, yr) {
    if (yr < min(tree_years) | yr > max(tree_years)) stop("invalid year")
    if (!(spc %in% sp_list$sp_name)) stop("invalid species")
    sp_attr <- filter(sp_list, sp_name == spc)
    min_ba <- sp_attr$rba_cm2 # basal area threshold
    subtree <- filter(trees, sp == str_to_lower(sp_attr$sp_code6))
    
    if (yr %in% tree_years) {
        subtree <- filter(subtree, census_year2 == yr, status2 == "A",
                          ba2 >= min_ba) %>%
            select(id = treeID, x, y, size = ba2)
    } else {  # interpolate growth and survival between census years
        subtree <- filter(subtree, status1 == "A", 
                           census_year1 < yr, census_year2 > yr) %>%
            mutate(ba_inter = ba_interp(ba1, ba2, census_year1, census_year2, yr),
                   st_inter = surv_interp(status2, census_year1, census_year2, yr)) %>%
            filter(ba_inter >= min_ba, st_inter == 1) %>%
            select(id = treeID, x, y, size = ba_inter)
    }
    subtree
}


# Subset "seedlings" data frame for a given species and year
#   buffer is the minimum distance from edge of plot
subset_seedlings <- function(spc, yr, buffer = 20) {
    if (!(yr %in% seedling_years)) stop("invalid year")
    if (!(spc %in% sp_list$sp_name)) stop("invalid species")
    # Only keep 'traps' that are not within buffer
    sub_traps <- filter(seedlings, PX > xmin + buffer, PX < xmax - buffer,
                        PY > ymin + buffer, PY < ymax - buffer)
    # Filter by species and remove status "P" (not in census yet)
    subseed <- filter(sub_traps, status != "P",
                      SPP %in% sp_list$sp_code6[sp_list$sp_name == spc])
    # Only keep records where seedling was first observed in given year
    subseed <- group_by(subseed, TAGF, PX, PY, SPP) %>%
        filter(year == min(year)) %>%
        ungroup() %>%
        filter(status == "A", year == yr)
    trap_counts <- group_by(subseed, PX, PY) %>%
        summarize(scount = n())
    # Need to add zeros for traps without seedlings
    trap_counts <- left_join(distinct(sub_traps, PX, PY), trap_counts)
    trap_counts$scount[is.na(trap_counts$scount)] <- 0
    rename(trap_counts, x = PX, y = PY, seeds = scount)
}


# Subset "seeds" data frame for species and year
#   buffer is the minimum distance from edge of plot
#   sel_traps is the subset of traps to select (or NULL if all traps)
subset_seeds <- function(spc, yr, buffer = 20, sel_traps = NULL) {
    if (!(spc %in% sp_list$sp_name)) stop("invalid species")
    # Subset by traps and excluding buffer
    sub_traps <- filter(seeds, X > xmin + buffer, X < xmax - buffer,
                               Y > ymin + buffer, Y < ymax - buffer)
    if (!is.null(sel_traps)) sub_traps <- filter(sub_traps, trap %in% sel_traps)
    # Subset by year and species
    subseed <- filter(sub_traps, year == yr,
                      sp %in% sp_list$sp_code4[sp_list$sp_name == spc])
    select(subseed, trap, x = X, y = Y, seeds = seed_eq)
}

    
# Subset "recruits" data frame for a given species and year
#   buffer is the minimum distance from edge of plot
#   sel_traps is the subset of (seed) traps to select (or NULL if all traps)
subset_recruits <- function(spc, yr, buffer = 20, sel_traps = NULL) {
    if (!(spc %in% sp_list$sp_name)) stop("invalid species")
    # Subset by traps and excluding buffer
    sub_traps <- distinct(select(seeds, trap, X, Y)) %>%
        filter(X > xmin + buffer, X < xmax - buffer,
               Y > ymin + buffer, Y < ymax - buffer)
    if (!is.null(sel_traps)) sub_traps <- filter(sub_traps, trap %in% sel_traps)
    sub_quads <- filter(recruits, trap %in% sub_traps$trap)
    # Subset by year (match seed year) and species
    subrec <- filter(sub_quads, seedyr == yr, 
                     sp %in% sp_list$sp_code4[sp_list$sp_name == spc])
    select(subrec, trap = quadrat, x = X, y = Y, seeds = recruits)
}

