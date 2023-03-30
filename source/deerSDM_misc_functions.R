sortBase <- function(vec, nsplines = 2, method, user_cp_quantiles = NULL, vecname) {
      #' Thin-plate spline basis function, modified from code
      #' provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
      #' penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1?24(2005).
      #'
      #' This function computes a thin-plate spline basis for a given numeric vector
      #' using a specified number of basis functions (nsplines).
      #' The basis functions are centered at control points determined by either the quantiles
      #' or kmeans clusters of the unique values in the input vector.
      #'
      #' @param vec A numeric vector for which the thin-plate spline basis is computed
      #' @param nsplines Integer specifying the number of basis functions (default: 2)
      #' @param methods Character string specifying the control point placement method
      #' @param user_cp_quantiles Optional numeric vector specifying user-defined quantiles for placing control points
      #' @return A matrix with the same number of rows as the input vector and nsplines columns,
      #'         representing the thin-plate spline basis for the input vector.

      if (!is.null(user_cp_quantiles)){
            control_points <- quantile(unique(vec), probs = user_cp_quantiles, na.rm = TRUE)
            plot_title <- paste0(vecname, ", Cp quantiles provided by user")
      } else {
            if (method == "quantile"){
                  control_points <- quantile(unique(vec), seq(0, 1, length = (nsplines+2))[-c(1, (nsplines+2))], na.rm = TRUE)
                  label = "quantiles"
            } else if (method == "kmeans"){
                  clusters <- kmeans(unique(na.omit(vec)), centers = nsplines)
                  control_points <- sort(clusters$centers)
                  label = "kmeans"
            } else if (method == "gmm"){
            require(mixtools)
            fit <- normalmixEM(unique(na.omit(vec)), k = nsplines)
            control_points <- fit$mu
            label = "Gaussian mixture model"
      }
            plot_title <- paste0(vecname, ", Cp based on ", label)
      }

      plot(density(na.omit(vec)), main = plot_title)
      abline(v = control_points, col = "red")

      # Define the thin-plate spline penalty matrix
      zFE       <- cbind(rep(1, length(vec) ), vec)   

      z_K <- (abs(outer(vec, control_points, "-"))) ^ 3
      OMEGA.all <- (abs(outer(control_points, control_points, "-"))) ^ 3
      svd.OMEGA.all  <- svd(OMEGA.all)
      sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u) *
                                                     sqrt(svd.OMEGA.all$d)))
      z.out     <- t(solve(sqrt.OMEGA.all, t(z_K)))
      return(z.out)
}

prepareSplineData <- function(df, vector, nsplines = 2, method = "quantile", user_cp_quantiles = NULL){
      #' Prepare spline data for a given vector in a data frame
      #'
      #' This function computes the spline basis for a given vector in a data frame using
      #' either natural splines or the custom sortBase method. The resulting spline basis
      #' is added to the data frame and prediction vectors are assigned to the global environment.
      #'
      #' @param df A data frame containing the input vector
      #' @param vector The vector within the data frame for which the spline basis is computed
      #' @param nsplines Integer specifying the number of basis functions (default: 2)
      #' @param methods Character string specifying the control point placement method (default: "quantile")
      #' @param user_cp_quantiles Optional numeric vector specifying user-defined quantiles for placing control points
      #' @return A modified data frame with added columns for each of the computed spline basis functions
      
      # Create prediction vector, as sequence of 100 values, spanning full range of values. 
      sequence <- seq(min(vector, na.rm = T), max(vector, na.rm = T), length = 100)
      # Save this vector, we use it later to label the x axis in the model effect plots
      assign(paste0(strsplit(deparse(substitute(vector)), split = "\\$")[[1]][2], ".seq"), sequence, envir = .GlobalEnv)

      # Add prediction vector to beginning of original values, scale all to mean 0 and sd 1
      var.s <- c(scale(c(sequence, vector)))  
      
      n_data <- NROW(vector)
      vec_name <- sub("^.*\\$", "", deparse(substitute(vector)))

      z.var.s <- scale(sortBase(vec = var.s, nsplines = nsplines, method, user_cp_quantiles, vecname = vec_name))
      for (i in 1:nsplines) {
            assign(paste0(strsplit(deparse(substitute(vector)), split = "\\$")[[1]][2], "_", i, ".sortBase"), c(z.var.s[1:100, i]), envir = .GlobalEnv)
            df[, paste0(strsplit(deparse(substitute(vector)), split = "\\$")[[1]][2], "_", i, ".sortBase")] <- c(z.var.s[101:(100+n_data), i])
      }

      
      return(df)
}

get_nearest_value <- function(main_vect, nearest_vect, value, max_dist = NULL, add_unique_value_name = NULL){
      #' Find the nearest values from a secondary vector for each element in the main vector
      #'
      #' This function finds the nearest value in the `nearest_vect` for each element in the
      #' `main_vect` based on their spatial location. Additional constraints and output
      #' modifications can be applied using the optional parameters.
      #'
      #' @param main_vect A spatial object representing the primary locations
      #' @param nearest_vect A spatial object representing the secondary locations
      #' @param value The name of the attribute in the nearest_vect to be assigned to the corresponding main_vect elements
      #' @param max_dist Numeric value specifying the maximum distance to consider for nearest points (default: NULL, no constraint)
      #' @param add_unique_value_name String to be appended to the attribute name in the output (default: NULL, no change)
      #' @return A data frame with the same number of rows as the main_vect and columns for x, y, and the nearest value(s) from nearest_vect
      if(is.lines(nearest_vect)) {
            nearest_vect <- as.points(nearest_vect, multi = TRUE)
      }
      closest_df <- as.data.frame(nearest(main_vect, nearest_vect))
      
      if(!is.null(max_dist)){
            closest_df <- closest_df %>% 
                  filter(distance <= max_dist)
      }
      
      na_df <- as.data.frame(main_vect[closest_df$from_id, ])
      na_df$x <- closest_df$from_x
      na_df$y <- closest_df$from_y
      na_df <- na_df[, c("x", "y")]
      if(!is.null(add_unique_value_name)){
            na_df[, paste0(value, "_new")] <- as.data.frame(nearest_vect[closest_df$to_id, ])[, value]
      }
      if(is.null(add_unique_value_name)){
            na_df[, value] <- as.data.frame(nearest_vect[closest_df$to_id, ])[, value]
      }
      
      na_df <- unique(na_df)
      return(na_df)
}

minute_difference <- function(tend, tstart) {
      tend <- as.character(tend)
      tstart <- as.character(tstart)
      
      tend <- sprintf("%04d", as.numeric(tend))
      tstart <- sprintf("%04d", as.numeric(tstart))
      
      tend <- paste0(substr(tend, 1, 2), ":", substr(tend, 3, 4))
      tstart <- paste0(substr(tstart, 1, 2), ":", substr(tstart, 3, 4))
      
      date_time_tend <- ymd_hm(paste0("1970-01-01 ", tend))
      date_time_tstart <- ymd_hm(paste0("1970-01-01 ", tstart))
      
      abs(as.numeric(difftime(date_time_tend, date_time_tstart, units = "mins")))
}

calc_vif <- function(data) {
      require(car)
      vif_values <- vif(lm(data))
      return(vif_values)
}



