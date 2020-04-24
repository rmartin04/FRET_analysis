## Libraries required
library(tidyverse)
library(here)

## Function required
## x is the input list, d is the distance allowed between coordinates between each time point
dfret <- function(x, d) {
  # Select columns of interest from dataframe
  y <- x %>%
    select('WellName', 'Field', 'Timepoint', 'X', 'Y', 'TrackedCells...FRET1', 'TrackedCells...CFP.Mean',
           'TrackedCells...YFP.Mean', 'TrackedCells...Cell.Roundness', 'TrackedCells...Cell.Area..µm..') %>%
  # Simplify column names
    rename('FRET1' = 'TrackedCells...FRET1',
           'CFP.Mean' = 'TrackedCells...CFP.Mean',
           'YFP.Mean' = 'TrackedCells...YFP.Mean', 
           'Cell.Roundness' = 'TrackedCells...Cell.Roundness',
           'Cell.Area.µm' = 'TrackedCells...Cell.Area..µm..')
  # Split dataframe by timepoint, now have list of dataframes, one for each time point
  y <- split(y, f = y$Timepoint) 
  
  # Loop through list of dataframes
  for (i in seq_along(y)){
    # Start with assigning a new object for the first time point
    if (i == 1){
      first <- y[[i]]
    }
    #For each time point, full_join by Field to match all objects between each time point
    else{
      first <- full_join(first, y[[i]], by = 'Field', suffix = c('', '.y')) %>%
        # Calculate the distance between the two points
        mutate(dist = sqrt((X.y - X)^2 + (Y.y - Y)^2)) %>%
        # Filter based off distance parameter, might be better to have multiple, relaxed parameters
        # than one strigent parameter
        filter(dist < d) %>%
        # Calculate dfret and df_fret
        mutate(dfret = FRET1.y - FRET1,
               df_fret = dfret / FRET1) %>%
        # Drop columns
        select(-c('X', 'Y', 'WellName.y', 'Timepoint.y')) %>%
        # Rename the new X, Y coordinates
        # Always going to be measuring distance between two consecutive time points
        # For each loop, the X and Y are updated
        rename('X' = 'X.y',
               'Y' = 'Y.y') %>%
        # Rename columns that were added with full_join to have the suffix as the time point
        rename_at(vars(ends_with('.y')), ~str_replace(., 'y', glue(i))) 
      # Rename other columns to add the time point suffix
      # Might be a better way to do this but I couldn't get it to work in the pipe
      names(first)[names(first) == 'dfret'] <- glue('dfret_', i)
      names(first)[names(first) == 'df_fret'] <- glue('df_fret_', i)
      names(first)[names(first) == 'dist'] <- glue('dist_', i)
    }
  }
  return(first)
}

# List all files in directory
files <- list.files(path = here(), pattern = "*].txt", full.names = T, recursive = T)

# Remove non-single cell files
akar_files <- files[grep("Population", files)]

# Read in files, skip if there is an error
akar_list <- lapply(akar_files, read.delim)

# Iterate through all files and apply the dfret function
# Alter the distance parameter
akar_fret <- lapply(akar, dfret, d = 50)

# Collapse list into one dataframe for further use
akar_fret <- plyr::ldply(akar_fret) 


