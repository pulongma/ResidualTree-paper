# Score all methods as Gaussian distributions with their predictive mean/variance.
# Full, PP and WN also supply their mixture for the separate *_exact scores.
# ResTree is required here even for the Vecchia and BART model files.
score_predictions <- function(predictions, prediction_object = NULL, y_center = 0) {
  gaussian_scores <- ResTree::restree_score(list(mean = predictions$mean, var = predictions$sd^2), predictions$observed,
    label = "gaussian"
  )
  scores <- data.frame(
    RMSE = gaussian_scores$RMSE,
    latent_RMSE = sqrt(mean((predictions$truth - predictions$mean)^2)),
    MAE = gaussian_scores$MAE,
    CRPS_Gaussian = gaussian_scores$CRPS,
    NLPD_Gaussian = gaussian_scores$logscore,
    coverage90_Gaussian = gaussian_scores$cover90,
    coverage95_Gaussian = gaussian_scores$cover95,
    width95_Gaussian = gaussian_scores$width95,
    z_var_Gaussian = gaussian_scores$z_var,
    CRPS_exact = NA_real_, NLPD_exact = NA_real_,
    coverage95_exact = NA_real_, width95_exact = NA_real_,
    score_type_exact = NA_character_
  )
  if (!is.null(prediction_object)) {
    # The tree models were fitted to the centered response.
    mixture_scores <- ResTree::restree_score(prediction_object, predictions$observed - y_center,
      label = "exact", exact = TRUE
    )
    scores$CRPS_exact <- mixture_scores$CRPS
    scores$NLPD_exact <- mixture_scores$logscore
    scores$coverage95_exact <- mixture_scores$cover95
    scores$width95_exact <- mixture_scores$width95
    scores$score_type_exact <- as.character(mixture_scores$score_type)
  }
  scores
}
