# Title: R source script for Beach geomorphology project
# Author: Michael Traurig

# Functions for Gravity analysis #

# Move points to nearest ocean cell. This will be used for both beach and port points
move_points_to_ocean <- function(points_sf, bathymetry, sitedata = NULL) {
  # Convert bathymetry to points, keeping only ocean cells
  for(i in 1:nrow(points_sf)){
    bath_croped <- crop(bathymetry, st_buffer(points_sf[i, ], 500))
    ocean_cells <- terra::as.points(bath_croped, values = TRUE, na.rm = TRUE)
    ocean_cells_sf <- st_as_sf(ocean_cells,
                               coords = c("x", "y"),
                               crs = st_crs(points_sf))
    if (!is.null(sitedata)) {
      colnames(ocean_cells_sf)[1] <- sitedata
    }
    if (i == 1) {
      points_sf_moved <- ocean_cells_sf[st_nearest_feature(points_sf[i, ],
                                                           ocean_cells_sf), ]
      if (!is.null(sitedata)) {
        points_sf_moved$Site[i] <- points_sf$Site[i]
      }
    } else {
      points_sf_moved <- rbind(points_sf_moved,
                               ocean_cells_sf[st_nearest_feature(points_sf[i, ],
                                                                 ocean_cells_sf), ])
      if (!is.null(sitedata)) {
        points_sf_moved$Site[i] <- points_sf$Site[i]
      }
    }
  }
  return(points_sf_moved)
}

## Find the nearest port to each study beach
find_nearest_port = function(beaches, ports) {

  get_nearest_portname = function(origin, destinations) {

    destinations[which.min(st_distance(origin, destinations)), ]$Site
  }

# apply to each beach, get_nearest_portname and assign the result to a new column in beaches
  for (i in 1:nrow(beaches)) {
    beaches$nearest_port[i] <- get_nearest_portname(beaches[i, ], ports)
  }
  return(beaches)
}

# Ocean transition function, converts bathyemtry and returns a transition layer
get_ocean_transition = function(bathy_raster, 	 #A bathymetry object
                                threshold = 3) { #Minimum depth threshold
  bathy_raster[bathy_raster >= threshold] <- NA
  bathy_raster <- abs(bathy_raster)
  #if bathy_raster is not a raster object, convert it to one
  if(!is(bathy_raster, "raster")) {
    bathy_raster <- raster(bathy_raster)
  }

	transOcean <- gdistance::transition(bathy_raster, mean, 8)
	gdistance::geoCorrection(transOcean, type="c")
}

# Ocean travel function
# This function takes two spatial points and a transition layer and returns the distance between the two points
get_ocean_travel = function(origin, #Start point
                            destination, #End point
                            transition_layer,
                            output = "distance") { #a transition layer object

  #if ocean_points is an sf object convert it to a spatial object
  if(is(origin, "sf")) {
    origin <- as(origin, "Spatial")
  }
destination <- filter(destination, Site == origin$nearest_port)
  if(is(destination, "sf")) {
    destination <- as(destination, "Spatial")
  }
  path <- st_as_sf(gdistance::shortestPath(transition_layer,
                                            origin,
                                            destination,
                                            output = "SpatialLines"))
  if(output == "distance") {
    return(st_length(path))
  } else {
    return(path)
  }

}
# Find nearest site function
# This function takes a spatial object of site locations and a single point spatial object and returns the nearest site to the point
get_nearest_site = function(site_points,		#a spatial object of site locations
							origin_point) {		#a single point spatial object for which the nearest site is to be found
  distances <- st_distance(origin_point, site_points)
  nearest_site_idx <- which.min(distances)
  nearest_site <- dplyr::select(sites[nearest_site_idx, ], "geometry")
  rbind(point_of_interest, nearest_site)
							}

get_land_distance = function(population, origin, structures, decay_factor = 0.0001, land_mask) {
  pop_reclass <- terra::classify(population, rcl = matrix(c(NA, NA, 0), ncol = 3, byrow = TRUE), right = NA) %>%
    mask(land_mask, touches = TRUE)

  distances <- distance(pop_reclass, origin) %>%
    mask(land_mask, touches = TRUE)

  road_raster <- terra::rasterize(vect(structures), pop_reclass, fun = "min", background = NA, touches = TRUE)
  distance_to_roads <- terra::distance(road_raster) %>%
    mask(land_mask, touches = TRUE)
  #Reclass distance to roads any values below 1 as 1
  dist_road <- terra::classify(distance_to_roads, rcl = matrix(c(0, 1, 1), ncol = 3, byrow = TRUE), right = NA)
  lambda <- decay_factor
  #decay_factor <- exp(-lambda * distance_to_roads)
  #distances_adj <- distances * decay_factor
  #decay_factor <- exp(lambda * (minmax(distance_to_roads)[2] - distance_to_roads))
  decay_factor <- exp(lambda * (minmax(dist_road)[2] - dist_road))
  distances_adj <- distances / decay_factor
  return(distances_adj)
}

calculate_gravity = function(study_site, population, constant, land_distances){
  total_distance <- land_distances[[which(names(land_distances) == study_site$nearest_port)]] + study_site$ocean_trav
  interactions <- constant * (values(population) * (1/values(total_distance)))
  gravity <- sum(interactions, na.rm = TRUE)
  return(gravity)
}
# Wave exposure functions
st_demtopoly = function(f,y) {
                  rast(f) %>% classify(matrix(c(minmax(.)[1]-1, y, NA, y, minmax(.)[2]+1, 1), ncol=3, byrow=TRUE)) %>% #bring in the raster and reclassify it
                  as.polygons() %>% st_as_sf() %>% st_make_valid() %>% dplyr::select()
                 }


# Debris data functions
Fix_2022_sites = function(dataset) {dataset %>% mutate(Event.ID = case_when(Site == Sites[1] ~ 1,
								   Site == Sites[2] ~ 2,
								   Site == Sites[3] ~ 3,
								   Site == Sites[4] ~ 4,
								   Site == Sites[5] ~ 5,
								   Site == Sites[6] ~ 6,
								   Site == Sites[7] ~ 7,
								   Site == Sites[8] ~ 8,
								   Site == Sites[9] ~ 9,
								   Site == Sites[10] ~ 10,
								   Site == Sites[11] ~ 11,
								   Site == Sites[12] ~ 12,
								   Site == Sites[13] ~ 13,
								   Site == Sites[14] ~ 14,
								   Site == Sites[15] ~ 15,
								   Site == Sites[16] ~ 16,
								   Site == Sites[17] ~ 17,
								   Site == Sites[18] ~ 18,
								   Site == Sites[19] ~ 19,
								   Site == Sites[20] ~ 20,
								   Site == Sites[21] ~ 21,
								   Site == Sites[22] ~ 22,
								   Site == Sites[23] ~ 23,
								   Site == Sites[24] ~ 24,
								   Site == Sites[25] ~ 25,
								   Site == Sites[26] ~ 26,
								   Site == Sites[27] ~ 27,
								   Site == Sites[28] ~ 28,
								   Site == Sites[29] ~ 29,
								   Site == Sites[30] ~ 30,
								   Site == Sites[31] ~ 31,
								   Site == Sites[32] ~ 32,
								   Site == Sites[33] ~ 33,
								   Site == Sites[34] ~ 34,
								   Site == Sites[35] ~ 35,
								   Site == Sites[36] ~ 36,
								   Site == Sites[37] ~ 37,
								   Site == Sites[38] ~ 38,
								   Site == Sites[39] ~ 39,
								   Site == Sites[40] ~ 40,
								   Site == Sites[41] ~ 41,
								   Site == Sites[42] ~ 42,
								   Site == Sites[43] ~ 43,
								   Site == Sites[44] ~ 44,
								   Site == Sites[45] ~ 45,
								   Site == Sites[46] ~ 46)) %>%
					mutate(Site = case_when(Site == "Saba.Bay..Hook.Island" ~ "Saba Bay, Hook Island",
								 Site == "Saba.Bay..Hook.Island.1" ~ "Saba Bay, Hook Island",
								 Site == "Saba.Bay..Hook.Island.2" ~ "Saba Bay, Hook Island",
								 Site == "Saba.Bay..Hook.Island.3" ~ "Saba Bay, Hook Island",
								 Site == "Haslewood.Island..Southern.End" ~ "Haselwood Island, Southern End",
								 Site == "haslewood.Island..Southern.End" ~ "Haselwood Island, Southern End",
								 Site == "Hazelwood.Island..Southern.End" ~ "Haselwood Island, Southern End",
								 Site == "Hazelwood.Island..Sotuhern.End" ~ "Haselwood Island, Southern End",
								 Site == "Hazelwood.Island..Sotuhern.End.1" ~ "Haselwood Island, Southern End",
								 Site == "Hazelwood.Island..Sotuhern.End.2" ~ "Haselwood Island, Southern End",
								 Site == "Mackerel.Bay.Hook.Island" ~ "Mackerel Bay Hook Island",
								 Site == "Mackerel.Bay.Hook.Island.1" ~ "Mackerel Bay Hook Island",
								 Site == "Mackerel.Bay.Hook.Island.2" ~ "Mackerel Bay Hook Island",
								 Site == "Mackeral.Bay..Hook.Island" ~ "Mackerel Bay Hook Island",
								 Site == "Cow.Reef..Cape.Conway.NP" ~ "Cow Reef, Cape Conway NP",
								 Site == "Long.Island.Sound" ~ "Long Island Sound",
								 Site == "Long.Island.Sound.1" ~ "Long Island Sound",
								 Site == "Maher.Island..Could.be.Maher.Island.East." ~ "Not sure",
								 Site == "Grimston.Point.Zone" ~ "Not sure",
								 Site == "Grimston.Point.Zone.1" ~ "Not sure",
								 Site == "Grimston.point.Zone" ~ "Not sure",
								 Site == "Turtle.Bay..Whitsunday.Island" ~ "Turtle Bay, Whitsunday Island",
								 Site == "Turtle.Bay..Whitsunday.Island.1" ~ "Turtle Bay, Whitsunday Island",
								 Site == "Turtle.Bay..Whitsunday.Island.2" ~ "Turtle Bay, Whitsunday Island",
								 Site == "Turtle.Bay..South.Molle.Island" ~ "Turtle Bay, South Molle Island",
								 Site == "Turtle.Bay..South.Molle.Island.1" ~ "Turtle Bay, South Molle Island",
								 Site == "Pandanus.Bay.Long.Island" ~ "Pandanus Bay Long Island",
								 Site == "Georges.Point" ~ "Georges Point",
								 Site == "Roberta.Bay..Shaw.Island" ~ "Roberta Bay, Shaw Island",
								 Site == "Roberta.Bay..Shaw.Island.1" ~ "Roberta Bay, Shaw Island",
								 Site == "Multiple.places.visited..cannot.use" ~ "Not sure",
								 Site == "Gloucester.Island..Eastern.Side" ~ "Gloucester Island, Eastern Side",
								 Site == "Airlie.Beach" ~ "Airlie Beach",
								 Site == "Multiple.locations.and.Not.specified.Maher.or.Maher.east" ~ "Not sure",
								 Site == "Cannot.use..two.locations" ~ "Not sure",
								 Site == "Cannot.use..two.locations.1" ~ "Not sure",
								 Site == "Cannot.use..two.locations.2" ~ "Not sure",
								 Site == "Cannot.use..Two.locations" ~ "Not sure",
								 Site == "Cannot.use..Two.locations.1" ~ "Not sure",
								 Site == "Cape.Conway" ~ "Cape Conway",
								 Site == "Cape.Conway.1" ~ "Cape Conway",
								 Site == "Cape.Conway.2" ~ "Cape Conway",
								 Site == "Cape.Conway.3" ~ "Cape Conway",
								 Site == "Cannot.use..Multiple.locations" ~ "Not sure",
								 Site == "Solace.Bay..Whitsunday.Island" ~ "Solace Bay, Whitsunday Island",
								 Site == "Thomas.Island..North.East.Beach" ~ "Thomas Island, North East Beach",
								 Site == "Unknown.location" ~ "Not sure")) %>%
								 mutate(year = 2022)
}


plot_overall = function(model, variable, xlab, ylab, title, limits, logged = TRUE) {
	datatoplot <- effects::effect(variable, model) %>% as.data.frame()
if(logged == TRUE) {
	ggplot(datatoplot, aes(x = get(variable), y = log(fit))) +
		geom_line() +
		geom_ribbon(aes(ymin = log(fit - se), ymax = log(fit + se)), alpha = 0.2) +
		labs(title = title,
			 x = xlab,
			 y = ylab) +
		ylim(limits) +
		theme_minimal() +
		theme(axis.title.y = element_text(size = 14),
		      axis.title.x = element_text(size = 14),
		  	  axis.text.x = element_text(size = 12),
         	  axis.text.y = element_text(size = 12),)
	} else {
	ggplot(datatoplot, aes(x = get(variable), y = fit)) +
		geom_line() +
		geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.2) +
		labs(title = title,
			 x = xlab,
			 y = ylab) +
		ylim(limits) +
		theme_minimal() +
		theme(axis.title.y = element_text(size = 14),
		      axis.title.x = element_text(size = 14),
		  	  axis.text.x = element_text(size = 12),
	     	  axis.text.y = element_text(size = 12),)
	}
}
plot_interaction = function(model, variable, xlab, title, ylab, limits, leg_title = 14, logged = TRUE, floatiness = "Unencumbered drift") {
 if(paste0("floatiness:",variable) %in% rownames(car::Anova(model)) == TRUE) {
	datatoplot <- effects::effect(paste0("floatiness:",variable), model) %>% as.data.frame()
  } else {
  datatoplot <- effects::effect(paste0(variable,":floatiness"), model) %>% as.data.frame()
  }
 if(floatiness == "Unencumbered drift") {
  } else {
	  datatoplot <- datatoplot %>% filter(!floatiness == "Unencumbered Drift")
  }
  if(logged == TRUE) {
	ggplot(datatoplot, aes(x = get(variable), y = log(fit), color = floatiness)) +
	  geom_line() +
	  geom_ribbon(aes(ymin = log(fit - se), ymax = log(fit + se)), alpha = 0.15, linetype = "blank") +
#	  facet_wrap(~floatiness) +
	  labs(title = title,
		   x = xlab,
		   y = ylab) +
	  ylim(limits) +
	  theme_minimal() +
 	  theme(legend.title = element_text(size = leg_title),
	  legend.text = element_text(size = 14),
	  axis.title.y = element_text(size = 14),
	  axis.title.x = element_text(size = 14),
  	  axis.text.x = element_text(size = 12),
	  axis.text.y = element_text(size = 12),
	  plot.title.position = "plot",
	  plot.title = element_text(hjust = 0.5))
  } else {
	ggplot(datatoplot, aes(x = get(variable), y = fit, color = floatiness)) +
	  geom_line() +
	  geom_ribbon(aes(ymin = fit - se, ymax = fit + se), alpha = 0.15, linetype = "blank") +
#	  facet_wrap(~floatiness) +
	  labs(title = title,
		   x = xlab,
		   y = ylab) +
	  ylim(limits) +
	  theme_minimal() +
	  theme(legend.title = element_text(size = leg_title),
	  legend.text = element_text(size = 14),
	  axis.title.y = element_text(size = 14),
	  axis.title.x = element_text(size = 14),
	  axis.text.x = element_text(size = 12),
	  axis.text.y = element_text(size = 12),
	  plot.title.position = "plot",
	  plot.title = element_text(hjust = 0.5))
  }
}

# ============================================================
# Function for predicting total debris across five categories
# ============================================================

five_floatiness_levels <- c(
  "Does Not Disperse",
  "Limited Dispersion",
  "Encumbered Drift",
  "Very Mobile",
  "Unencumbered Drift"
)

predict_five_category_debris <- function(
    new_data,
    fitted_model,
    scaling_parameters) {

  # Remove geometry temporarily if an sf object was supplied.

  model_inputs <- if (
    inherits(new_data, "sf")
  ) {
    st_drop_geometry(new_data)
  } else {
    as.data.frame(new_data)
  }

  required_columns <- c(
    "Site",
    "gravity",
    "site_exposure"
  )

  missing_columns <- setdiff(
    required_columns,
    names(model_inputs)
  )

  if (length(missing_columns) > 0) {
    stop(
      paste(
        "The following prediction variables are missing:",
        paste(
          missing_columns,
          collapse = ", "
        )
      )
    )
  }

  model_inputs <- model_inputs %>%
    dplyr::select(
      Site,
      gravity,
      site_exposure
    ) %>%
    mutate(
      Site = as.character(Site),
      gravity = as.numeric(gravity),
      site_exposure = as.numeric(site_exposure)
    )

  if (anyDuplicated(model_inputs$Site)) {
    stop(
      paste(
        "Each row supplied to the prediction function must",
        "have a unique Site identifier."
      )
    )
  }

  if (
    anyNA(model_inputs) ||
    any(!is.finite(model_inputs$gravity)) ||
    any(!is.finite(model_inputs$site_exposure))
  ) {
    stop(
      "Missing or non-finite predictor values were supplied."
    )
  }

  # Generate five rows for each coastline location.

  prediction_data <- tidyr::crossing(
    model_inputs,
    floatiness = five_floatiness_levels
  ) %>%
    mutate(
      floatiness = factor(
        floatiness,
        levels = five_floatiness_levels
      ),
      num_cleans = 1
    )

  # A valid training-year level is supplied because year appears
  # in the model structure. It does not affect these predictions
  # because re.form = NA excludes the year random effect.

  reference_year <- as.character(
    fitted_model$frame$year[1]
  )

  training_year_levels <- levels(
    fitted_model$frame$year
  )

  if (is.null(training_year_levels)) {
    training_year_levels <- unique(
      as.character(
        fitted_model$frame$year
      )
    )
  }

  prediction_data <- prediction_data %>%
    mutate(
      Site = factor(Site),
      year = factor(
        reference_year,
        levels = training_year_levels
      )
    )

  # Generate population-level predictions.
  #
  # re.form = NA excludes both Site and year random effects.
  # Predictions therefore represent the relationships with
  # gravity, exposure and debris dispersal category.

  prediction_data$predicted_debris <-
    predict(
      fitted_model,
      newdata = prediction_data,
      type = "response",
      re.form = NA,
      allow.new.levels = TRUE
    )

  # Retain category-specific predictions.

  category_predictions <-
    prediction_data %>%
    dplyr::select(
      Site,
      floatiness,
      predicted_debris
    )

  # Sum the five predictions to obtain one accumulation score
  # for each location.

  total_predictions <-
    category_predictions %>%
    group_by(Site) %>%
    summarise(
      debris = sum(
        predicted_debris
      ),
      .groups = "drop"
    )

  list(
    totals = total_predictions,
    categories = category_predictions
  )
}
