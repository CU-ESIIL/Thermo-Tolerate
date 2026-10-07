# Download ERA5 for the bounding box of all Rothamsted trap sites, one request
# per year, in parallel. Safe to stop and re-run: months whose .nc already
# exists are skipped, and a failed year doesn't stop the others.
# Run from the 3_theory project (e.g. RStudio Background Job, or
# `Rscript scripts/download_era5.R` from 3_theory/).

library(sf)
library(mcera5)
library(here)
library(future)
library(furrr)

out_path <- here("data", "mcera5")
years    <- 1990:2023

roth_coords <- read.csv(here("data", "rothamsted_coordinates.csv"))
bbox <- st_bbox(st_multipoint(cbind(roth_coords$longitude, roth_coords$latitude)))
bbox <- bbox + c(-0.25, -0.25, 0.25, 0.25)  # one ERA5 cell of padding for edge/coastal sites

download_year <- function(y) {
  # build_era5_request() returns one sub-request per month
  req <- build_era5_request(xmin = bbox[["xmin"]], xmax = bbox[["xmax"]],
                            ymin = bbox[["ymin"]], ymax = bbox[["ymax"]],
                            start_time = as.POSIXct(paste0(y, "-01-01 00:00"), tz = "UTC"),
                            end_time   = as.POSIXct(paste0(y, "-12-31 23:00"), tz = "UTC"),
                            outfile_name = paste0("era5_roth_", y))

  # A month is done once its unzipped .nc exists (a leftover .zip may be partial)
  month_nc <- file.path(out_path, sub("\\.zip$", ".nc", sapply(req, `[[`, "target")))
  todo     <- req[!file.exists(month_nc)]

  if (length(todo) > 0) {
    # overwrite = TRUE only replaces leftover .zip files for months still to do;
    # request_era5() otherwise errors on any existing file instead of skipping it
    res <- tryCatch(
      request_era5(request = todo, out_path = out_path, overwrite = TRUE, combine = FALSE),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      message("Year ", y, " failed: ", conditionMessage(res), " -- re-run to resume")
      return(invisible(FALSE))
    }
  }

  # Combine the 12 months into one file per year, matching mcera5's own naming
  year_nc <- file.path(out_path, paste0("era5_roth_", y, "_", y, ".nc"))
  if (all(file.exists(month_nc)) && !file.exists(year_nc)) {
    combine_netcdf(filenames = month_nc, combined_name = year_nc)
  }
  invisible(all(file.exists(month_nc)))
}

# CDS only runs a few requests per user at once, so more workers just queue
plan(multisession, workers = 4)
done <- future_map_lgl(years, download_year,
                       .options = furrr_options(seed = TRUE, scheduling = Inf,
                                                packages = c("sf", "mcera5")))
plan(sequential)

message("Complete years: ", sum(done), " of ", length(years))
if (any(!done)) message("Still missing: ", paste(years[!done], collapse = ", "))
