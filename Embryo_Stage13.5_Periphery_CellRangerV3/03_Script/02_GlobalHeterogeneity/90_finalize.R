# ########################################################
# This script wraps up analysis (clean, save, render, ...)
# ########################################################


## @knitr final_saveSessionImage

# Save an image of the current R session containing all objects
sessionImagePath = file.path( PATH_OUTPUT, paste0( outputFilesPrefix, "sessionImage.RDATA"));
save(list = ls( all.names = TRUE), 
     file = sessionImagePath, 
     envir = environment());

#save.image( sessionImagePath); # Saves everything in .Globalenv which is not what we want depending on rendering context



