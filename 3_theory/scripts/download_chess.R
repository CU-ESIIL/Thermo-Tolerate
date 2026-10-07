# Download CHESS-met v1.2 (1961-2019, 1 km daily, Great Britain) from the EIDC
# and keep only the grid cell at each Rothamsted trap site.
#
# EIDC has no server-side subsetting, so each full-GB monthly file is downloaded,
# the trap cells are cut out, and the full file is deleted. Disk use stays tiny;
# download volume does not.
#
# Output, one folder per site, ready for micro_uk(spatial = <site folder>):
#   data/micro_uk/sites/<Trap>/chess_<var>_<YYYYMM>.nc   (1 x 1 cell)
#   data/micro_uk/sites/<Trap>/terr50.tif, terr1000.tif  (cropped from data/micro_uk)
# One cell per site also avoids micro_uk() choosing the wrong cell: it matches
# longitude along the grid's southern edge, which can be ~8 km off in the east.
#
# Before running:
#   1. Accept the CHESS-met licence on the EIDC catalogue page
#      https://catalogue.ceh.ac.uk/id/835a50df-e74f-4bfb-b593-804fd61d5eab
#   2. Create a personal access token (https://eidc.ac.uk/help/personalaccesstokens)
#      and add EIDC_TOKEN=<token> to your user .Renviron (usethis::edit_r_environ()),
#      never to a file in this repository. Restart R.
#   3. Build data/micro_uk/terr50.tif and terr1000.tif (OS Terrain 50).
# Safe to stop and re-run: finished months are skipped.

library(here)
library(terra)
library(ncdf4)
library(future)
library(furrr)

years    <- 1990:2019
vars     <- c("dtr", "tas", "huss", "precip", "rsds", "sfcWind", "psurf")  # what micro_uk() reads
dem_dir  <- here("data", "micro_uk")
site_dir <- here("data", "micro_uk", "sites")
dataset  <- "https://catalogue.ceh.ac.uk/datastore/eidchub/835a50df-e74f-4bfb-b593-804fd61d5eab"

# File layout on the EIDC datastore. If the test download below fails with 404,
# check one file's link on the EIDC download page and edit this function.
chess_url <- function(var, y, m) {
  last_day <- format(seq(as.Date(sprintf("%d-%02d-01", y, m)), by = "month", length.out = 2)[2] - 1, "%d")
  sprintf("%s/%s/chess-met_%s_gb_1km_daily_%d%02d01-%d%02d%s.nc", dataset, var, var, y, m, y, m, last_day)
}

token <- Sys.getenv("EIDC_TOKEN")
if (token == "") stop("Set EIDC_TOKEN in your user .Renviron (see header) and restart R")

download_chess <- function(url, dest) {
  h <- curl::new_handle()
  curl::handle_setheaders(h, Authorization = paste("Bearer", token))
  curl::curl_download(url, dest, handle = h, quiet = TRUE)
}

# --- Sites and their CHESS cells --------------------------------------------
sites <- read.csv(here("data", "rothamsted_coordinates.csv"))
sites <- sites[!duplicated(sites[c("longitude", "latitude")]), ]  # Kirton / Kirton 2 share a point
sites$id <- gsub("[^A-Za-z0-9]+", "_", trimws(sites$Trap))
pts <- project(vect(sites, geom = c("longitude", "latitude"), crs = "EPSG:4326"), "EPSG:27700")
sites[c("easting", "northing")] <- crds(pts)

# Use one month of tas to fix each site's cell; coastal sites whose own cell is
# sea (NA) move to the nearest land cell
test_file <- file.path(tempdir(), "chess_test.nc")
message("Test download: ", chess_url("tas", years[1], 1))
download_chess(chess_url("tas", years[1], 1), test_file)
nc <- nc_open(test_file)
gx <- ncvar_get(nc, "x"); gy <- ncvar_get(nc, "y")
tas1 <- ncvar_get(nc, "tas", start = c(1, 1, 1), count = c(-1, -1, 1))  # x by y, first day
nc_close(nc)
land <- which(!is.na(tas1), arr.ind = TRUE)
for (s in seq_len(nrow(sites))) {
  d <- (gx[land[, 1]] - sites$easting[s])^2 + (gy[land[, 2]] - sites$northing[s])^2
  k <- which.min(d)
  sites$ix[s] <- land[k, 1]; sites$iy[s] <- land[k, 2]
  sites$moved_m[s] <- round(sqrt(d[k]))  # distance from site to its cell centre
}
print(sites[c("Trap", "easting", "northing", "moved_m")])  # > ~700 m means the site's own cell is sea

# --- Cropped elevation rasters per site -------------------------------------
terr50   <- rast(file.path(dem_dir, "terr50.tif"))
terr1000 <- rast(file.path(dem_dir, "terr1000.tif"))
for (s in seq_len(nrow(sites))) {
  d <- file.path(site_dir, sites$id[s]); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  e <- ext(sites$easting[s] - 3000, sites$easting[s] + 3000, sites$northing[s] - 3000, sites$northing[s] + 3000)
  if (!file.exists(file.path(d, "terr50.tif")))   writeRaster(crop(terr50, e),   file.path(d, "terr50.tif"))
  if (!file.exists(file.path(d, "terr1000.tif"))) writeRaster(crop(terr1000, e), file.path(d, "terr1000.tif"))
}
write.csv(sites, file.path(site_dir, "site_cells.csv"), row.names = FALSE)

# --- Is a per-site file complete? ------------------------------------------
# A file only counts as done if it opens, holds the right variable as a
# 1 x 1 x <days in month> array, and is not all NA. Existence alone is not
# enough: interrupted runs left missing files and full-GB files under final names.
days_in_month <- function(y, m) {
  as.integer(format(seq(as.Date(sprintf("%d-%02d-01", y, m)), by = "month", length.out = 2)[2] - 1, "%d"))
}
site_file_status <- function(f, var, y, m) {
  if (!file.exists(f)) return("missing")
  if (file.size(f) > 1e5) return("full_grid")   # a 1-cell month is ~1 KB
  nc <- tryCatch(ncdf4::nc_open(f), error = function(e) NULL)
  if (is.null(nc)) return("unreadable")
  on.exit(ncdf4::nc_close(nc))
  if (!var %in% names(nc$var)) return("wrong_var")
  if (!identical(as.integer(nc$var[[var]]$varsize), c(1L, 1L, days_in_month(y, m)))) return("wrong_dims")
  if (all(is.na(ncdf4::ncvar_get(nc, var)))) return("all_NA")
  "ok"
}
site_file_ok <- function(f, var, y, m) site_file_status(f, var, y, m) == "ok"

# --- Cut one site's cell out of a full-GB file into micro_uk's file name ----
write_site_cell <- function(full, var, s, out) {
  nc   <- nc_open(full)
  vals <- ncvar_get(nc, var, start = c(s$ix, s$iy, 1), count = c(1, 1, -1))
  tvals <- ncvar_get(nc, "time")
  tatt <- ncatt_get(nc, "time")
  vatt <- ncatt_get(nc, var)
  nc_close(nc)

  dx <- ncdim_def("x", "m", gx[s$ix])
  dy <- ncdim_def("y", "m", gy[s$iy])
  dt <- ncdim_def("time", tatt$units, tvals, unlim = TRUE, calendar = if (is.null(tatt$calendar)) "standard" else tatt$calendar)
  v  <- ncvar_def(var, if (is.null(vatt$units)) "" else vatt$units, list(dx, dy, dt), missval = 1e20)
  tmp <- paste0(out, ".part")
  o <- nc_create(tmp, v)
  ncvar_put(o, v, array(vals, c(1, 1, length(vals))))
  nc_close(o)
  file.rename(tmp, out)  # only complete files get the final name, so reruns can trust it
}

# --- Download every month in parallel ---------------------------------------
months <- expand.grid(m = 1:12, y = years)

process_month <- function(y, m) {
  for (var in vars) {
    outs <- file.path(site_dir, sites$id, sprintf("chess_%s_%d%02d.nc", var, y, m))
    good <- vapply(outs, site_file_ok, logical(1), var = var, y = y, m = m)
    if (all(good)) next
    unlink(outs[!good])  # remove bad/partial files so they are rebuilt
    full <- file.path(tempdir(), sprintf("chess_%s_%d%02d_full.nc", var, y, m))
    ok <- tryCatch({
      download_chess(chess_url(var, y, m), full)
      ncdf4::nc_close(ncdf4::nc_open(full))  # fail here if the download is truncated
      TRUE
    }, error = function(e) { message(var, " ", y, "-", m, " download failed: ", conditionMessage(e)); FALSE })
    if (ok) for (s in which(!good)) {
      tryCatch(write_site_cell(full, var, sites[s, ], outs[s]),
               error = function(e) message(var, " ", y, "-", m, " ", sites$id[s], " cut failed: ", conditionMessage(e)))
    }
    unlink(full)
  }
  # Done only if every variable is valid at every site
  all(vapply(vars, function(var) {
    all(vapply(file.path(site_dir, sites$id, sprintf("chess_%s_%d%02d.nc", var, y, m)),
               site_file_ok, logical(1), var = var, y = y, m = m))
  }, logical(1)))
}

plan(multisession, workers = 4)
done <- future_map2_lgl(months$y, months$m, process_month,
                        .options = furrr_options(seed = TRUE, scheduling = Inf,
                                                 packages = c("ncdf4", "curl")))
plan(sequential)

message("Complete months: ", sum(done), " of ", nrow(months), " -- re-run to retry any failures")
