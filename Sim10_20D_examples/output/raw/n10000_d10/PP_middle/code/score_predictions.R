# Score all methods as Gaussian distributions with their predictive mean/variance.
# Full, PP and WN also supply their mixture for the separate *_exact scores.
# ResTree is required here even for the Vecchia and BART model files.
score_predictions <- function(predictions, prediction_object = NULL, y_center = 0) {
  p <- predictions
  g <- ResTree::restree_score(list(mean = p$mean, var = p$sd^2), p$observed,
                              label = "gaussian")
  scores <- data.frame(
    RMSE = g$RMSE,
    latent_RMSE = sqrt(mean((p$truth - p$mean)^2)),
    MAE = g$MAE,
    CRPS_Gaussian = g$CRPS,
    NLPD_Gaussian = g$logscore,
    coverage90_Gaussian = g$cover90,
    coverage95_Gaussian = g$cover95,
    width95_Gaussian = g$width95,
    z_var_Gaussian = g$z_var,
    CRPS_exact = NA_real_, NLPD_exact = NA_real_,
    coverage95_exact = NA_real_, width95_exact = NA_real_,
    score_type_exact = NA_character_)
  if (!is.null(prediction_object)) {
    # The tree models were fitted to the centered response.
    e <- ResTree::restree_score(prediction_object, p$observed - y_center,
                                label = "exact", exact = TRUE)
    scores$CRPS_exact <- e$CRPS
    scores$NLPD_exact <- e$logscore
    scores$coverage95_exact <- e$cover95
    scores$width95_exact <- e$width95
    scores$score_type_exact <- as.character(e$score_type)
  }
  scores
}
