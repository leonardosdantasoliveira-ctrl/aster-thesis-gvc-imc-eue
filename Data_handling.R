library(tidyverse)
library(data.table)
library(countrycode)

## 1. Main MRIO and GVC variables (Load Data) ----

start.year = "2000"
end.year = "2019"

start_mrio <- read.csv(
  paste(start.year,".csv",sep=""),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

end_mrio <- read.csv(
  paste(end.year,".csv",sep=""),
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mrios <- list(start_mrio,end_mrio)  # get mrios of start and end year in a list with the years as the names, so now every time we do an lapply or Map function it will apply everything inside the function equally to the start year and end year data
names(mrios) <- c(start.year,end.year)

ghg <- fread("E.csv",
             select = c("Emissions origin area", "ACTIVITY", "TIME_PERIOD", "OBS_VALUE"))
setnames(ghg, c("country", "parent", "year", "value")) # called sector parent for now as will be useful in manipulation below
ghg$country <- ifelse(ghg$country == "Rest of the world","ROW",countrycode(ghg$country, origin = "country.name", destination = "iso3c")) # import environmental extension data which we will format in greater detail after fixing the MRIOs


# first function removes non Z and Y rows and colummns and separates the Z and Y tables
mrios <- lapply(mrios, function(df) {
  MRIO <- as.matrix(df)
  MRIO <- MRIO[-which(row.names(df) %in% c("TLS","VA","OUT")),-which(colnames(df) %in% "OUT")] # removes non Z and Y rows
  mode(MRIO) <- "numeric"
  
  n_cs_local <- nrow(MRIO)                                    
  Z <- MRIO[, 1:n_cs_local, drop = FALSE]                     # Z is square: n_cs × n_cs
  Y <- MRIO[, (n_cs_local + 1):ncol(MRIO), drop = FALSE]      # Y is whatever columns come after
  
  list(
    Z = Z, # saves consolidated Z matrix in mrios list
    Y = Y # saves Y matrix in mrios list
  )
})

ghg_icio_mismatch <- unique(sub("_.*$", "", row.names(mrios[[start.year]]$Z))[which(!sub("_.*$", "", row.names(mrios[[start.year]]$Z)) %in% ghg$country)]) # gets all the countries that are in the Z matrix but not present in the environmentally extended data

mrios <- lapply(mrios, function(df) {  # function that will combine our extra MRIO countries into the ROW in the MRIO data
  Z <- df[["Z"]] # locally defines our Z table from the mrios list so we can use Z in this function and not df[["Z"]] every time
  Y <- df[["Y"]] # locally defines our Y  table from the mrios list so we can use Y in this function and not df[["Y"]] every time
  row_secs <- sub("^ROW_", "", rownames(Z)[grepl("^ROW_", rownames(Z))]) # specific sectors of ROW's sectors 
  row_fds <- sub("^ROW_", "", colnames(Y)[grepl("^ROW_", colnames(Y))]) # specific column titles of ROW's final demand
  for (cty in ghg_icio_mismatch) { # loop that goes, sector by sector, and final demand component by final demand component, through each of our extra countries' values across the Z and Y tables and cumulatively adds them to ROW's values
    cty_labs <- paste0(cty, "_", row_secs) # lookup in the Z table for the country in the loop that is to be combined with ROW's sectors
    row_labs <- paste0("ROW_", row_secs) # lookup in the Z table for ROW's sectors
    cty_fd_cols <- paste0(cty, "_", row_fds) # lookup in the Y table for the country in the loop that is to be combined with ROW's final demand components
    row_fd_cols <- paste0("ROW_", row_fds) # lookup in the Y table for ROW's final demand components
    Z[row_labs, ] <- Z[row_labs, ] + Z[cty_labs, ] # combines all the given country's rows of the Z table of the same sector with all ROW's rows of the Z table, in other words the rows corresponding to ROW's sectors are increased one country at a time in the loop
    Y[row_labs, ] <- Y[row_labs, ] + Y[cty_labs, ] # combines all the given country's rows of the Y table of the same sector with all ROW's rows of the Y table, in other words the rows corresponding to ROW's sectors are increased one country at a time in the loop
    Z[, row_labs] <- Z[, row_labs] + Z[, cty_labs] # combines all the given country's columns of the Z table of the same sector with all ROW's columns of the Z table, in other words the columns corresponding to ROW's sectors are increased one country at a time in the loop
    Y[, row_fd_cols] <- Y[, row_fd_cols] + Y[, cty_fd_cols] # combines all the given country's columns of the Y table of the same final demand component with all ROW's columns of the Y table, in other words the columns corresponding to ROW's sectors are increased one country at a time in the loop
  }
  drop_pattern <- paste0("^(", paste(ghg_icio_mismatch, collapse = "|"), ")_") # defines the country prefixes to drop from the dataset now that these mismatched countries are fully integrated into ROW
  Z <- Z[!grepl(drop_pattern, rownames(Z)), !grepl(drop_pattern, colnames(Z))] # keeps all rows and columns of the Z matrix that don't have the mismatched country prefixes, in other words drops our combined into ROW countries from the Z matrix
  Y <- Y[!grepl(drop_pattern, rownames(Y)), !grepl(drop_pattern, colnames(Y))] # keeps all rows and columns of the Y matrix that don't have the mismatched country prefixes, in other words drops our combined into ROW countries from the Y matrix
  df$Z <- Z # save our newly condensed Z table in our mrios list, overwriting the old one
  df$Y <- Y # save our newly condensed Y table in our mrios list, overwriting the old one
  df  # important that this is here so the new version df is returned to the mrios list that includes what we just added to it; if this line isn't here the new Z and Y tables won't return to the mrios list
})

## 1.1Set sizes ----
n_cs <- nrow(mrios[[start.year]][["Z"]])  # total number of country-sectors in the dataset (should obviously be the same if it's start.year or end.year)
row_country <- substr(row.names(mrios[[end.year]][["Z"]]), 1, 3) # the full country order of the dataset
country_codes <- unique(row_country) # specific countries in the data set
n_countries_z <- length(country_codes) # number of countries in the data set
n_sec <- n_cs / n_countries_z # number of sectors per country in the dataset
n_fd_per_country <- ncol(mrios[[start.year]][["Y"]]) / n_countries_z # number of final demand components per country in the dataset

row_labs     <- rownames(mrios[[start.year]][["Z"]]) # full country-sector order of the data set, should be the same whether it's start.year or end.year
col_labs     <- colnames(mrios[[start.year]][["Z"]]) # full country-sector order of the data set, should be the same whether it's start.year or end.year, indeed should be the same as row_labs as the Z matrix is symmetrical in country-sector order which(row_labs!=col_labs) should return 0
fd_col_labs  <- colnames(mrios[[start.year]][["Y"]]) # full country/final demand component order of the Y matrix
fd_countries <- substr(fd_col_labs, 1, 3) # full country order of the Y matrix, useful when matching a particular country with its 6 columns in the Y matrix

sector_code <- sub("^[^_]+_", "", col_labs) # removes country codes from labels so we are only left with sectors and we can identify primary product sectors or digital / intangible sectors along the columns
target_sectors <- unique(sector_code) # list of full 50 sectors present in the data set, target sectors to keep when we filter imported labor share or environmental extension data
is_primary <- grepl("^(A|B)", sector_code) | sector_code == "C19" # A or B as primary product sectors, plus C19 (refined petroleum) treated as primary due to proximity to extraction
digital_sectors <- c("J61", "J62_63") # definition of digital sectors
intangible_sectors <- c("J58T60", "J61", "J62_63", "M") # definition of intangible sectors 
imc_core_set <- c("USA", "GBR", "DEU", "FRA", "NLD", "CHE", "IRL", "JPN", "KOR", "TWN", "CHN") # definition of 'core' IMC countries
is_digital_row    <- sector_code %in% digital_sectors 
is_intangible_row <- sector_code %in% intangible_sectors 

mrios <- lapply(mrios, function(df) { # function we need to get the global Leontief inverse and vectors of output, value-added, and value-added per output that will be useful throughout the rest of the script
  Z <- df[["Z"]] # locally defines our Z table from the mrios list so we can use Z in this function and not df[["Z"]] every time
  Y <- df[["Y"]] # locally defines our Y table from the mrios list so we can use Y in this function and not df[["Y"]] every time
  x <- rowSums(Z) + rowSums(Y) # total output per country-sector (should be the same as the OUT column in the originally imported file we deleted in the first function)
  A <- matrix( # sets up dimensions of the A matrix to fill in below
    0,
    nrow = n_cs,
    ncol = n_cs,
    dimnames = list(rownames(Z), colnames(Z))
  )
  positive_x <- x > 0 # to avoid dividing by 0 when output is 0
  A[, positive_x] <- sweep(Z[, positive_x, drop = FALSE], 2, x[positive_x], "/")  #  this is the technical coefficients matrix, A
  B_global <- solve(diag(n_cs) - A) # this is the global Leontief inverse, (I - A)^-1
  va <- x - colSums(Z) # VA by industry = output - intermediate inputs (could also take the VA row we previously deleted from the originally imported data but sometimes there's a minor amount of taxes in a tax row that gets subtracted from the total output and slightly reduces VA and that's avoided this way)
  v_coeff <- rep(0, n_cs) # sets the dimensions of the value-added / output vector to be filled in below
  v_coeff[positive_x] <- va[positive_x] / x[positive_x]  # gets VA / output vector
  v_coeff[abs(v_coeff) < 1e-12] <- 0  # makes extremely small values 0
  list( # saves all our new objects in the mrios list, important that Z and Y are also here or else they won't be saved even though we didn't change them inside this function
    Z = Z,
    Y = Y,
    x = x,
    va = va,
    A = A,
    B_global = B_global,
    v_coeff = v_coeff
  )
})

labr_path <- "LABR.csv"  

labr <- fread(labr_path,
              select = c("REF_AREA", "ACTIVITY", "COUNTERPART_AREA",
                         "MEASURE",  "UNIT_MEASURE", "TIME_PERIOD", "OBS_VALUE"))
setnames(labr, c("country", "sector", "counterpart",
                 "measure", "unit", "year", "value"))

labr <- labr[counterpart == "W"         &
               measure     == "LABR"      &
               unit        == "PT_VA"]  # import labor share data, taking the compensation of employees as a share of VA measure
labr$sector <- dplyr::recode(labr$sector,"C241_2431" = "C24A","C242_2432" = "C24B") # so that sectors in the labor share data are compatible name-wise with the MRIO sectors
labr <- labr[sector %in% target_sectors] # keep the level of sector disaggregation that is the same as what we have with the MRIO data, i.e. 50 sectors
labr <- labr[country %in% country_codes] # keep the countries that exist in our MRIO dataset, get rid of aggregates like WXOECD, etc

labr <- labr %>% # due to some outlier values, we will keep 3 years surrounding the start and end year in order to replace implausible values with nearby year averages
  mutate(key = paste(country, sector, sep = "_")) %>%
  filter(year %in% c(as.numeric(start.year):(as.numeric(start.year) + 3),
                     (as.numeric(end.year) - 3):as.numeric(end.year))) %>%
  mutate(
    year_col = paste0("value_", year)
  ) %>%
  select(key, country, sector, year_col, value) %>%
  pivot_wider(names_from = year_col, values_from = value) %>%
  mutate(
    "{paste0('value_', start.year)}" := if_else(
      .data[[paste0("value_", start.year)]] < 10 |
        .data[[paste0("value_", start.year)]] > if_else(sector == "T", 100, 90),
      rowMeans(pick(all_of(paste0("value_", as.numeric(start.year) + 1:3))), na.rm = TRUE),
      .data[[paste0("value_", start.year)]]
    ),
    "{paste0('value_', end.year)}" := if_else(
      .data[[paste0("value_", end.year)]] < 10 |
        .data[[paste0("value_", end.year)]] > if_else(sector == "T", 100, 90),
      rowMeans(pick(all_of(paste0("value_", as.numeric(end.year) - 3:1))), na.rm = TRUE),
      .data[[paste0("value_", end.year)]]
    )
  ) %>%
  bind_rows( # adds extra rows from ROW which is not in the labor share data set so that the length of the data frame is the same as the length of all the MRIO data
    expand.grid(
      country = country_codes[!country_codes %in% labr$country],
      sector = target_sectors,
      stringsAsFactors = FALSE
    ) %>%
      mutate(key = paste(country, sector, sep = "_"))
  ) %>%
  arrange(match(key, names(mrios[[start.year]]$va))) # assure that the country-sector order is the same as in the MRIO data set

split_map <- list( # to split sectors in the environmentally extended data where the dataset has higher aggregate categories than the MRIO data
  A01_02    = c("A01", "A02"),
  B05_06    = c("B05", "B06"),
  B07_08    = c("B07", "B08"),
  C24       = c("C24A", "C24B"),
  C30       = c("C301", "C302T309")
)

ghg <- ghg %>% filter(year==start.year|year==end.year) %>% filter(parent %in% target_sectors|parent %in% names(split_map)) %>% filter(country %in% country_codes) %>% mutate(year = as.character(year)) # some preliminary cleaning of the GHG dataset, keeping the right level of country and sector level detail to match the MRIO data as closely as possible but still needs to split a few aggregated sectors

mrios <- Map(function(df, yr) { # this function further treats the environmentally extended and labor share data to put it in its final form, yr is defined as names(mrios) i.e. the start year and the end year
  x <- df[["x"]] # locally defines our output vector from the mrios list so we can use x in this function and not df[["x"]] every time
  va <- df[["va"]] # locally defines our value-added vector from the mrios list so we can use x in this function and not df[["x"]] every time
  ghg_yr <- ghg %>% dplyr::filter(year == yr) # gets environmentally extended data in year specific format
  direct <- ghg_yr[!parent %in% names(split_map),
                   .(country, icio_sector = parent, value)] # keeps the data for all sectors that match 1:1 the MRIO data and therefore don't need to be split
  x_dt <- data.frame("row_labs" = names(x), "country" = substr(names(x), 1, 3), "sector" = sub("^[^_]*_", "", names(x)), "x" = x) # gets a country-sector data frame of total output from the MRIOs for all country-sectors to be used below to split aggregated environmentally extended sectors
  x_dts <- list() # to store results of the loop below looping through each aggregated sector from the environmentally extended data
  for(i in 1:length(split_map)){ # loop through each aggregated sector from the environmentally extended data to split emissions by share of output represented by each MRIO disaggregated sector in the environmentally extended aggregated category
    x_dts[[i]] <- x_dt %>% filter(sector %in% split_map[[i]]) %>% mutate(parent = names(split_map)[i]) %>% group_by(country, parent) %>% mutate(total = sum(x, na.rm = TRUE),share = if_else(total > 0, x / total, 0)) %>% ungroup() %>% inner_join(.,ghg_yr,by = c("country","parent")) %>% mutate(value = value * share) %>% dplyr::rename("icio_sector" = sector) %>% select(country,icio_sector,value)
  }
  splits <- dplyr::bind_rows(x_dts) # put the data stored separately in the loop above together so we now have all country-sector data on emissions that needed to be split
  ghg_long <- rbind(direct,splits) # add previously stored 1:1 sector matches with newly split data to get the full data set of emissions for all country-sectors in the dataset
  ghg_long <- ghg_long %>% mutate(key = paste(country,icio_sector,sep="_")) %>% arrange(key) # need to bring back output vector which we dropped to merge with direct data frame which didn't have an output vector, which is done immediately below
  ghg_long$x <- x[ghg_long$key] # merge output vector with full post-split ghg data in a way that is certain that row order isn't an issue
  df$ghg_long <- ghg_long %>% mutate(f = if_else(x > 0, value / x, 0)) # create the full emissions intensity vector and store it in the mrios list
  labr <- labr %>% select(key, country, sector, all_of(paste0("value_", yr))) # get the labor share data in a year specific way to be treated within the function, dropping the other years
  labr$va <- va[labr$key] # add a value-added column to the labor share data in a way that is certain that row order isn't an issue
  labr <- labr %>% dplyr::rename(coe_labshare = all_of(paste0("value_", yr))) %>% # rename the colname of the value we want coe_labshare
    mutate( # above outside of the function we replaced implausible values with surrounding 3 year averages, but since there are still some remaining implausible values, we now want to ultimately fix upper and lower plausible limits to all values
      coe_labshare = if_else(
        !is.na(coe_labshare),
        pmin(
          pmax(coe_labshare, 10),
          if_else(sector == "T", 100, 90)
        ),
        coe_labshare
      ),
      coe_labshare = if_else(
        !is.na(coe_labshare) & va == 0,
        0,
        coe_labshare
      )
    )
  df$coe_labshare <- labr # save final labor share data in the mrios list
  df # important that this is here so the new version df is returned to the mrios list that includes what we just added to it; if this line isn't here the new Z and Y tables won't return to the mrios list
}, mrios, names(mrios)) # defines df and yr as mrios and names(mrios) in the arguments of the function

results_list <- vector("list", length = length(country_codes)) # to store country-specific results in the below function
mrios <- lapply(mrios, function(df) { # one more function this time to get all remaining MRIO based variables for the dataset
  ghg_long <- df[["ghg_long"]] # locally defines ghg_long data frame from the mrios list so we can use ghg_long in this function and not df[["ghg_long"]] every time
  Z <- df[["Z"]] # locally defines our Z table from the mrios list so we can use Z in this function and not df[["Z"]] every time
  Y <- df[["Y"]] # locally defines our Y table from the mrios list so we can use Y in this function and not df[["Y"]] every time
  B_global <- df[["B_global"]] # locally defines our B_global table from the mrios list so we can use B_global in this function and not df[["B_global"]] every time
  labr <- df[["coe_labshare"]] # locally defines our coe_labshare data frame from the mrios list so we can use coe_labshare in this function and not df[["coe_labshare"]] every time
  Z_colsums <- colSums(Z) # sums all intermediate purchases across all country sectors domestic and foreign
  digital_int_by_col <- colSums(Z[is_digital_row, , drop = FALSE]) # sums all digital sector intangible purchases across all country sectors domestic and foreign
  strict_denom_digital <- colSums(Z[!is_digital_row, , drop = FALSE]) # sums all NON digital sector intangible purchases across all country sectors domestic and foreign
  intan_int_by_col <- colSums(Z[is_intangible_row, , drop = FALSE]) # sums all intangible sector intangible purchases across all country sectors domestic and foreign
  strict_denom_intan <- colSums(Z[!is_intangible_row, , drop = FALSE]) # sums all NON intangible sector intangible purchases across all country sectors domestic and foreign
  Intangibles_Dependency_dig       <- ifelse(Z_colsums > 0, digital_int_by_col / Z_colsums, NA_real_) # gets Intangibles_Dependency_dig var, share of total domestic and foreign intermediate purchases that are to digital sectors
  Intangibles_Dependency_dig_strict <- ifelse(Z_colsums > 0, digital_int_by_col / strict_denom_digital, NA_real_) # gets Intangibles_Dependency_dig_strict var, share of total domestic and foreign intermediate purchases that are to digital sectors (with digital sector purchases removed from denominator)
  Intangibles_Dependency_intan       <- ifelse(Z_colsums > 0, intan_int_by_col / Z_colsums, NA_real_) # gets Intangibles_Dependency_intan var, share of total domestic and foreign intermediate purchases that are to intangible sectors
  Intangibles_Dependency_intan_strict <- ifelse(Z_colsums > 0, intan_int_by_col / strict_denom_intan, NA_real_) # gets Intangibles_Dependency_intan_strict var, share of total domestic and foreign intermediate purchases that are to intangible sectors (with intangible sector purchases removed from denominator)
  FID_num_any   <- numeric(n_cs)   # empty vector to be filled later with foreign intangible inputs, any source
  FID_num_imc   <- numeric(n_cs)   # empty vector to be filled later foreign intangible inputs, from IMC core
  FID_num_usa   <- numeric(n_cs)   # empty vector to be filled later foreign intangible inputs, from USA only
  Z_foreign_any <- numeric(n_cs)   # empty vector to be filled later all foreign intermediate inputs (alt denom.)
  f <- ghg_long$f[match(row_labs, ghg_long$key)] # defines locally the emissions intensity vector, f, in a way that is certain to be in the same row order as the MRIO data
  country_index <- split(seq_along(row_labs),sub("_.*$", "", row_labs))[country_codes] # index of country rows in the MRIO data set row order
  Y_country <- matrix(0, nrow = n_cs, ncol = n_countries_z,  # set empty vector to fill with aggregated final demand per country
                      dimnames = list(row_labs, country_codes))
  for(cc in country_codes){ # loop through all countries and aggregate final demand by country 
    Y_country[, cc] <- rowSums(Y[, fd_countries == cc, drop = FALSE])
  }
  CO2_footprint_country <- as.numeric(matrix(f, nrow = 1) %*% B_global %*% Y_country) # country based material footprint measure, all emissions domestic and foreign embodied in a country's final demand purchases
  
  # (Not used anymore) Getting sector-based material footprint: total global emissions embodied in each country-sector's final-demand output
   
  fd_total              <- rowSums(Y)                              # n_cs vector
  CO2_footprint_sector  <- as.numeric(matrix(f, nrow = 1) %*% B_global) * fd_total
  names(CO2_footprint_sector) <- row_labs
  
  # Alternative sectoral standardizer: same f %*% B_global formula, but normalized
  # to total output (x) rather than final demand (Y). 
  x_total                 <- rowSums(Z) + rowSums(Y)
  CO2_total_chain_x       <- as.numeric(matrix(f, nrow = 1) %*% B_global) * x_total
  names(CO2_total_chain_x) <- row_labs
  
  names(CO2_footprint_country) <- colnames(Y_country)
  for(cc in country_codes){  # loop through all countries to get the variables that require country-specific local vectors or matrices, they will be stored in the country-specific results_list within the loop
    idx <- country_index[[cc]] # define country-specific rows in the MRIO row order
    B_cc <- B_global[idx, idx, drop = FALSE]  # extracts the local Leontief inverse of a particular country's sectors
    V_c <- matrix(df$v_coeff[idx], nrow = 1) # extracts V_c (the VA / output shares) for all sectors in a particular country
    Vc_hat <- diag(as.numeric(V_c)) # diagonalizes V_c
    Ehat_c <- diag(f[idx]) # gets the country specific emissions intensities for all sectors in a particular country and diagonalizes
    f_foreign      <- f   
    f_foreign[idx] <- 0 # 0s out the domestic rows for a particular country of the emissions intensity vector, thus creating f_foreign
    foreign_int_cols <- setdiff(seq_len(nrow(Z)), idx) # identifies rows in the Z table that are foreign to our particular country
    exgr_intermediate <- rowSums(Z[idx, foreign_int_cols, drop = FALSE]) # row sums for gross intermediate exports from Z for all of our particular country's sectors
    foreign_fd_cols <- which(fd_countries != cc) # identify all final demand columns that are not from our particular country
    exgr_final <- rowSums(Y[idx, foreign_fd_cols, drop = FALSE]) # row sums for gross final demand exports from Y for all of our particular country's sectors
    exgr_sector <- exgr_intermediate + exgr_final # gives us a vector of total gross exports for all of the sectors from our particular country, which is EXGRc
    exgr_full <- matrix(exgr_sector, ncol = 1)
    sector_code_local <- sub("^[^_]+_", "", names(exgr_sector))  
    is_primary_local <- grepl("^(A|B)", sector_code_local) | sector_code_local == "C19" # as we did above globally, but locally this time, including C19 as primary
    exgr_vec_np <- exgr_sector # duplicate EXGRc in order to set primary product sectors to 0
    exgr_vec_np[is_primary_local] <- 0 # set primary products to 0 for gross exports vector for our particular country
    exgr_vec_np <- matrix(exgr_vec_np, ncol = 1) # make sure our non primary product version of EXGRc vector is a matrix so it's ready to multiply with matrix multiplication, operation will be performed a bit further below
    nonprimary_rows <- row_labs[(!is_primary)] # removes primary product rows from all countries from the Z matrix 
    foreign_nonprimary_rows <- nonprimary_rows[which(!grepl(cc,nonprimary_rows))] # removes all remaining domestic industry rows from our particular country from the Z matrix, leaving us with the rows we want to column sum only non primary gross intermediate imports for all columns, which will be iPM_np, which we get a bit further below
    fn_idx <- which(!is_primary & substr(row_labs, 1, 3) != cc) # same as foreign_nonprimary_rows but in indexed numbers and not row label names, easier to use on empty matrices that may not have row labels
    M_full <- matrix(0, nrow = n_cs, ncol = length(idx)) # empty matrix to store Z colsums that are foreign to our particular country
    M_full[foreign_int_cols, ] <- Z[foreign_int_cols, idx, drop = FALSE] # storing them
    M_np <- matrix(0, nrow = n_cs, ncol = length(idx)) # empty matrix to store Z colsums that are foreign to our particular country and exclusive of primary products from all countries
    M_np[fn_idx, ] <- Z[fn_idx, idx, drop = FALSE] # storing them
    foreign_intang_rows <- which(is_intangible_row & row_country != cc) # index to identify intangible sector rows that are foreign to our particular country
    foreign_intang_imc  <- which(is_intangible_row & row_country %in% setdiff(imc_core_set, cc)) # index to identify intangible sector rows that are foreign to our particular country and from one of the earlier defined 'core' IMC countries
    foreign_intang_usa  <- which(is_intangible_row & row_country == "USA" & "USA" != cc)  # index to identify intangible sector rows that are foreign to our particular country and from the USA alone (will be empty when our particular country is the USA)
    foreign_all_rows    <- which(row_country != cc) # index to identify all foreign rows 
    FID_num_any[idx]   <- colSums(Z[foreign_intang_rows, idx, drop = FALSE]) # numerator for FID_any, sums total intangible intermediate imports for all of our particular country's sectors that are foreign to our particular country
    FID_num_imc[idx]   <- colSums(Z[foreign_intang_imc,  idx, drop = FALSE]) # numerator for FID_imc, sums total intangible intermediate imports for all of our particular country's sectors that are foreign to our particular country and from the 'core' IMC countries
    FID_num_usa[idx]   <- colSums(Z[foreign_intang_usa,  idx, drop = FALSE]) # numerator for FID_imc, sums total intangible intermediate imports for all of our particular country's sectors that are from the USA (unless our particular country is the USA)
    Z_foreign_any[idx] <- colSums(Z[foreign_all_rows,    idx, drop = FALSE]) # denominator for FID measures, total foreign intermediate imports for our particular country
    CO2_footprint_sector_c <- CO2_footprint_sector[idx]
    CO2_total_chain_x_c <- CO2_total_chain_x[idx]
    
    dva_total <- as.numeric(Vc_hat %*% B_cc %*% exgr_vec_np) # our formula to get XDVA_np (non primary products) for all sectors in our country, that is, total domestic value added accruing to each sector from total gross exports of any sector in the economy, except for primary product exports
    imint_np <- colSums(Z[foreign_nonprimary_rows, which(grepl(cc,col_labs)), drop = FALSE]) # gives us iPM_np, a vector of our non primary gross intermediate imports for all sectors in our particular country from the Z matrix 
    XDCO2_full <- as.numeric(Ehat_c %*% B_cc %*% exgr_full) # environmental version of XDVA except not exclusive of primary products
    XDCO2_np   <- as.numeric(Ehat_c %*% B_cc %*% exgr_vec_np) # environmental version of XDVA that is excluding primary products
    MCO2_full <- as.numeric(matrix(f_foreign, nrow = 1) %*% (B_global %*% M_full)) # environmental version of iPM except not exclusive of primary products
    MCO2_np   <- as.numeric(matrix(f_foreign, nrow = 1) %*% (B_global %*% M_np)) # environmental version of iPM except that is excluding primary products
    net_imp_np   <- MCO2_np   - XDCO2_np # environmental balance variable (non-standardized, level) for the primary products excluded version
    net_imp_full <- MCO2_full - XDCO2_full # environmental balance variable (non-standardized, level) for the version with primary products included
    net_pct_country_np <- 100 * net_imp_np / CO2_footprint_country[cc] # environmental balance variable standardized by country footprint
    net_pct_country_full <- 100 * net_imp_full / CO2_footprint_country[cc] # environmental balance variable standardized by country footprint
    net_pct_sector_np   <- ifelse(CO2_footprint_sector_c > 0,
                                  100 * net_imp_np   / CO2_footprint_sector_c,
                                  NA_real_) # environmental balance variable standardized by country-sector footprint
    net_pct_sector_full <- ifelse(CO2_footprint_sector_c > 0,
                                  100 * net_imp_full / CO2_footprint_sector_c,
                                  NA_real_) # environmental balance variable standardized by country-sector footprint
    net_pct_chain_x_np   <- ifelse(CO2_total_chain_x_c > 0,
                                   100 * net_imp_np   / CO2_total_chain_x_c,
                                   NA_real_)
    net_pct_chain_x_full <- ifelse(CO2_total_chain_x_c > 0,
                                   100 * net_imp_full / CO2_total_chain_x_c,
                                   NA_real_)
    results_list[[cc]] <- data.frame(
      country                 = cc,
      industry                = row_labs[idx],            
      va                      = df$va[idx],               # millions of USD
      XDVA_np = dva_total,                                # millions of USD
      iPM_np = imint_np,                                  # millions of USD
      E_direct                = f[idx],                    # millions of tonnes CO2e
      CO2_footprint_country   = CO2_footprint_country[cc],  # millions of tonnes CO2e (country-level, repeated)
      CO2_footprint_sector  = CO2_footprint_sector_c,
      CO2_total_chain_x         = CO2_total_chain_x_c,
      
      ## non-primary version 
      CO2_exp_np              = XDCO2_np,                      # millions of tonnes CO2e
      CO2_imp_np              = MCO2_np,                       # millions of tonnes CO2e
      CO2_net_imp_np          = net_imp_np,                    # millions of tonnes CO2e
      
      ## full version (for robustness, includes primary products)
      CO2_exp_full            = XDCO2_full,                # millions of tonnes CO2e
      CO2_imp_full            = MCO2_full,                # millions of tonnes CO2e
      CO2_net_imp_full        = net_imp_full,              # millions of tonnes CO2e
      CO2_net_imp_pct_chain_full = net_pct_chain_x_full,
     
      
      ## main indicator + robustness
      CO2_net_imp_pct_np      = net_pct_country_np,            # (% country footprint, non-primary)
      CO2_net_imp_pct_full    = net_pct_country_full,          # robustness (% country footprint, full)
      CO2_net_imp_pct_np_sec   = net_pct_sector_np,           # (% country-sector footprint, non-primary)
      CO2_net_imp_pct_full_sec = net_pct_sector_full,          # robustness (% country-sector footprint, full)
      CO2_net_imp_pct_chain_np   = net_pct_chain_x_np,        # (% country-sector footprint based on x instead of FD, non-primary)
      CO2_net_imp_pct_chain_full = net_pct_chain_x_full,       # robustness (% country-sector footprint, full)
      
      stringsAsFactors = FALSE
    )
  }
  
  FID_share_total <- ifelse(Z_colsums > 0, FID_num_any / Z_colsums, NA_real_) # share of foreign intangible inputs in total intermediate inputs
  FID_share_foreign <- ifelse(Z_foreign_any > 0, FID_num_any / Z_foreign_any, NA_real_) # share of foreign intangible inputs in total foreign intermediate inputs
  
  
  
  ##### Added to get the necessary inputs for the construction of IMC Peripherality
  # Getting foreign share of total intangible inputs
  foreign_intan_share_total_intan <- ifelse(intan_int_by_col > 0,
                                            FID_num_any / intan_int_by_col,
                                            NA_real_)
  
  # Deviation from global mean (negative = below-average foreign share = more 'core',
  #      positive = above-average foreign share = more 'periphery')
  foreign_intan_share_dev <- foreign_intan_share_total_intan -
    median(foreign_intan_share_total_intan, na.rm = TRUE)
  
  # IMC core–periphery index, weighted by importance of intangibles
  IMC_positioning <- foreign_intan_share_dev * Intangibles_Dependency_intan
  
  ## Rent-concentration within the foreign intangible basket:
  FID_imc_share_of_fid <- ifelse(FID_num_any > 0, FID_num_imc / FID_num_any, NA_real_) # after the loop we now have the full country-sector vectors of FID_any and FID_imc and can thus divide them to get FID_imc_share_of_fid
  FID_usa_share_of_fid <- ifelse(FID_num_any > 0, FID_num_usa / FID_num_any, NA_real_) # after the loop we now have the full country-sector vectors of FID_any and FID_usa and can thus divide them to get FID_usa_share_of_fid (note this will be 0 for USA sectors)
  #######
  
  
  
  imc_vars <- data.frame(
    row_country, sector_code, row_labs,
    Intangibles_Dependency_dig, Intangibles_Dependency_dig_strict, Intangibles_Dependency_intan, Intangibles_Dependency_intan_strict,
    FID_num_any, FID_num_imc, FID_num_usa,
    FID_share_total, FID_share_foreign,
    FID_imc_share_of_fid, FID_usa_share_of_fid,
    IMC_positioning
  )
  
  colnames(imc_vars) <- c(
    "country", "sector", "key",
    "Intangibles_Dependency_dig", "Intangibles_Dependency_dig_strict", "Intangibles_Dependency_intan", "Intangibles_Dependency_intan_strict",
    "FID_num_any", "FID_num_imc", "FID_num_usa",
    "FID_share_total", "FID_share_foreign",
    "FID_imc_share_of_fid", "FID_usa_share_of_fid",
    "IMC_positioning"
  )
  
  df$imc_vars <- imc_vars # store imc vars on their own in mrios
  df$results_list <- results_list #store individual country-based list of key variables on their own in mrios
  results <- dplyr::bind_rows(results_list) # put all individually stored variables together to get a full cs length data frame

  results$coe_labshare <- labr$coe_labshare[match(results$industry, labr$key)] # add labor share value to the full cs length data frame in a way that is certain to match the row order of the full cs length data frame
  results <- results %>%  
    left_join(imc_vars %>% select(-c(country)),
              by = c("industry" = "key")) %>%
    mutate(GVC_part = (XDVA_np + iPM_np) / va,
           GVC_capt = XDVA_np / (XDVA_np + iPM_np),
           CO2_net_imp_np_per_va = ifelse(va > 0,
                                          CO2_net_imp_np / va,
                                          NA_real_)) %>%
    group_by(country) %>%
    mutate(va_share = va / sum(va, na.rm = TRUE)) %>%
    ungroup()
  df$results <- results # save full dataset results in mrios as results
  df # important that this is here so the new version df is returned to the mrios list that includes what we just added to it within this function
})

dataset <- purrr::imap(mrios, function(df, yr) { # combines the start year and end year results to form one integrated dataset
  df$results %>%
    dplyr::rename_with(
      ~ paste0(.x, "_", yr), # appends _XXXX year to each variable depending on whether they came from start or end year results
      -c(country, sector, industry)
    )
}) %>%
  purrr::reduce(  # joins the two datasets 
    dplyr::full_join,
    by = c("country", "sector", "industry")
  ) %>%
  select(country,sector,industry,everything())

## Cleaning structurally-empty rows for analysis ----

# Diagnostic: how many rows have zero VA in both years (= sector doesn't exist for this country)?
empty_both <- with(dataset, va_2000 == 0 & va_2019 == 0)
cat("Rows with zero VA in both years:", sum(empty_both, na.rm = TRUE), "\n")

# Diagnostic: sector T is residual and consistently produces NAs
cat("Sector T rows:", sum(dataset$sector == "T"), "\n")

# Build the analysis-ready sample
# Excluded sectors:
#   T       — residual category (households as employers), near-zero economic content
#   A01-A03 — primary: agriculture, forestry, fishing — not GVC-relevant
#   B05-B09 — primary: mining and extraction — not GVC-relevant
#   C19     — refined petroleum: behaves as primary, close to extraction
#   L       — real estate: non-tradable, non-GVC

excluded_sectors <- c(
  "T",
  paste0("A0", 1:3),
  paste0("B0", 5:9),
  "C19",
  "L"
)

dataset_full     <- dataset                                  # keep the full version with NAs
dataset_analysis <- dataset %>%
  filter(!sector %in% excluded_sectors) %>%                  # drop residual + primary + L
  filter(!(va_2000 == 0 & va_2019 == 0))                     # drop sectors that don't exist anywhere

cat("Excluded sectors:     ", paste(excluded_sectors, collapse = ", "), "\n")
cat("Full dataset rows:    ", nrow(dataset_full), "\n")
cat("Analysis-ready rows:  ", nrow(dataset_analysis), "\n")

# Per-variable NA counts in the cleaned sample
na_summary <- dataset_analysis %>%
  summarise(across(
    c(starts_with("CO2_net_imp_pct"),
      starts_with("FID_"),
      starts_with("Intangibles_Dependency_")),
    ~ sum(is.na(.))
  )) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_NA") %>%
  arrange(desc(n_NA))

print(na_summary, n = Inf)

# Save both
write.csv(dataset_full,     "Data_2000_2019_full.csv",     row.names = FALSE)
write.csv(dataset_analysis, "Data_2000_2019_analysis.csv", row.names = FALSE)




