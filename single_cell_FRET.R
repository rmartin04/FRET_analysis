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
      full <- y[[i]]
    }
    #For each time point, full_join by Field to match all objects between each time point
    else{
      filter <- full_join(full, y[[i]], by = c('Field'), suffix = c('', '.y')) %>%
        # Calculate the distance between the two points
        mutate(dist.y = sqrt((X.y - X)^2 + (Y.y - Y)^2)) %>%
        # Filter based off distance parameter, might be better to have multiple, relaxed parameters
        # than one strigent parameter
        filter(dist.y < d) %>%
        # Group_by cell coordinates of preceding time point
        group_by(X, Y) %>%
        # select pairing with closest pairings
        top_n(-1, dist.y) %>%
        # Group_by cell coordinates of matched time point
        group_by(X.y, Y.y) %>%
        # select pairing with closest pairings, these top_n make sure that each cell is only kept once
        # in case it matches to two cells within the d parameter
        top_n(-1, dist.y) %>%
        # Calculate dfret and df_fret
        mutate(dfret.y = FRET1.y - FRET1,
               df.fret.y = dfret.y / FRET1)
      
      # Select columns from the matched time point
      matched <- filter[,str_detect(colnames(filter), ".y")]
      # Add back the coordinates from the previous time point
      matched$X <- filter$X
      matched$Y <- filter$Y
      
      # Left join by the coordinates from the previous time point
      # This will add NAs for cells that have not been matched
      full <- left_join(full, matched, by = c('X', 'Y'))
      # rename previous coordinates with the timepoint suffix
      names(full)[names(full) == 'X'] <- paste('X', i - 1, sep = '_')
      names(full)[names(full) == 'Y'] <- paste('Y', i - 1, sep = '_')
      # replace coordinates used for full_join in the next iteration
      names(full)[names(full) == 'X.y'] <- 'X'
      names(full)[names(full) == 'Y.y'] <- 'Y'
      
      # drop unnecessary columns
      full$Timepoint.y <- NULL
      full$WellName.y <- NULL
      
      # Replace '.y' suffix of matched columns with the timepoint
      full <- full %>%
        rename_at(vars(ends_with('.y')), ~str_replace(., '.y', paste0("_", i))) 
    }
  }
  # The last time point matched will not have the timepoint suffix added, this adds it
  len <- length(y)
  names(full)[names(full) == 'X'] <- paste('X', len, sep = "_")
  names(full)[names(full) == 'Y'] <- paste('Y', len, sep = "_")
  # Add the timepoint suffix to the first time point columns
  names(full)[names(full) == 'FRET1'] <- paste('FRET1', 1, sep = "_")
  names(full)[names(full) == 'CFP.Mean'] <- paste('CFP.Mean', 1, sep = "_")
  names(full)[names(full) == 'YFP.Mean'] <- paste('YFP.Mean', 1, sep = "_")
  names(full)[names(full) == 'Cell.Roundness'] <- paste('Cell.Roundness', 1, sep = "_")
  names(full)[names(full) == 'Cell.Area.µm'] <- paste('Cell.Area.µm', 1, sep = "_")
  
  # Pivot table to tidy format 
  full <- full %>%
    rownames_to_column() %>%
    pivot_longer(cols = 5:ncol(.)) %>%
    separate(name, into = c("axis", "timepoint"), sep = "_") %>% 
    drop_na() %>%
    pivot_wider(names_from = axis, values_from = value) %>%
    mutate(df.fret = replace_na(df.fret, 0),
           dfret = replace_na(dfret, 0))
  
  # Return the table
  return(full)
}

# List all files in directory
files <- list.files(path = here(), pattern = "*].txt", full.names = T, recursive = T)

# Remove non-single cell files
akar_files <- files[grep("Population", files)]

# Read in files
akar <- map(akar_files, read.delim)

# Iterate through all files and apply the dfret function
# Alter the distance parameter
akar_fret <- map_dfr(akar, dfret, d = 50)


