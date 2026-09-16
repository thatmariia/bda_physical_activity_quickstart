# Modelling functions

#' Get feature data excluding outcomes and confidence cols
#' @param df The input data frame
#' @return A data frame with only the feature columns
df_feat <- function(df) {
    df |> select(-aggr_activity, -activity_confidence)
}

#' Get a data frame with feature cols and additional
#' cols for distances between data points and k-means centers
#' @param df The input data frame
#' @param km The k-means object
#' @return A data frame with the k-means distances
get_kmeans_dist_df <- function(df, km) {
    dists <- proxy::dist(df_feat(df), km$centers)
    df_km_dists <- unclass(dists) |>
        data.frame() |>
        mutate(aggr_activity = df$aggr_activity)
    return(df_km_dists)
}

#' Get a data frame with feature cols and additional
#' col for the k-means cluster assignment
#' @param df The input data frame
#' @param km The k-means object
#' @return A data frame with the k-means cluster assignments
get_kmeans_df <- function(df, km) {
    df_km <- data.frame(
        cluster = factor(km$cluster),
        aggr_activity = df$aggr_activity
    )
    return(df_km)
}

#' Define options for model training, which are used to
#' construct different model configurations
#' @return A data frame with the training options
define_train_options <- function() {
    opts <- expand_grid(
        method = "multinom",
        term_mode = c("nomode", "inter", "add"),
        term = c("noterm", "km_cluster", "km_dist"),
        weighting = c("unweighted", "weighted")
    ) |>
        filter(!(term_mode == "nomode" & term != "noterm")) |>
        filter(!(term == "noterm" & term_mode != "nomode"))
    return(opts)
}

#' Construct a list of specifications for fitting a model
#' based on provided options
#' @param df The input data frame
#' @param df_km The data frame with k-means cluster assignments
#' @param df_km_dists The data frame with k-means distances
#' @param method The method for training
#' @param term_mode The term mode option
#' @param term The term option
#' @param weighting The weighting option
#' @return A list of specifications for fitting a model
get_fit_spec <- function(
    df, df_km, df_km_dists,
    method, term_mode, term, weighting
) {
    key <- paste(method, term_mode, term, weighting, sep = "_")

    data <- df_feat(df) |>
        mutate(aggr_activity = df$aggr_activity) |>
        mutate(
            weights = if (weighting == "weighted") {
                df$activity_confidence
            } else {
                rep(1, nrow(df))
            }
        )

    # add terms
    data <- switch(term,
        "noterm" = data,
        "km_cluster" = data |> bind_cols(cluster = df_km$cluster),
        "km_dist" = data |> bind_cols(df_km_dists |> select(-aggr_activity)),
        data
    )

    # add term mode
    if (term_mode == "inter" && term == "km_cluster") {
        formula <- aggr_activity ~ (. - weights) * cluster
    } else if (term_mode == "inter" && term == "km_dist") {
        feat_names <- names(df_feat(df))
        dist_names <- names(df_km_dists |> select(-aggr_activity))
        # ==> START LLM https://chatgpt.com/share/6aaab81f-4d5c-83eb-bfac-021950015cd3
        rhs <- paste0(
            "(", paste(feat_names, collapse = " + "), ") * ",
            "(", paste(dist_names, collapse = " + "), ")"
        )
        formula <- as.formula(paste("aggr_activity ~", rhs))
        # ==> END LLM
    } else {
        formula <- aggr_activity ~ (. - weights)
    }

    return(list(data = data, formula = formula, method = method, key = key))
}

#' Preprocess and fit a list of models based on the provided data and options
#' @param df The input data frame
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @param nstart The number of random starts for k-means
#' @return A list of fitted models
fit_models <- function(df, number = 2, repeats = 1, nstart = 2) {
    # Define train control with repeated cross-validation
    trcntr <- caret::trainControl(method = "repeatedcv", number = number, repeats = repeats, verboseIter = FALSE)

    # Precompute relevant data (with additional features)
    k <- length(unique(df$aggr_activity))
    km_fit <- kmeans(df_feat(df), centers = k, nstart = nstart)
    df_km <- get_kmeans_df(df, km_fit)
    df_km_dists <- get_kmeans_dist_df(df, km_fit)

    # Define options for training
    opts <- define_train_options()

    # Get specifications for training models
    specs <- pmap(opts, \(method, term_mode, term, weighting) {
            get_fit_spec(
                df, df_km, df_km_dists,
                method, term_mode, term, weighting
            )
        }
    )

    # Fit the models
    models <- list()
    for (spec in specs) {
        cat("Fitting model:", spec$key, "...\n")
        model <- caret::train(
            spec$formula,
            data = spec$data,
            weights = weights,
            method = spec$method,
            trControl = trcntr,
            preProcess = c("nzv", "corr", "center", "scale"),
            MaxNWts = 10000,
            maxit = 1000,
            trace = FALSE
        )
        models[[spec$key]] <- model
    }

    return(models)
}

#' Fit all models and return the results
#' @param df The input data frame
#' @param weights The weights for the model (NULL for default)
#' @param number The number of folds for cross-validation
#' @param repeats The number of times to repeat the cross-validation
#' @param nstart The number of random starts for k-means
#' @return A list of fitted models and their results
fit_all <- function(df, number = 2, repeats = 1, nstart = 2) {
    models <- fit_models(df, number, repeats, nstart)

    results <- data.frame(
        accuracy = sapply(models, function(x) max(x$results$Accuracy)),
        kappa = sapply(models, function(x) max(x$results$Kappa))
    )

    return(list(models = models, results = results))
}
