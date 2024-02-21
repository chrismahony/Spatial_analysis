#functions
coloc_one_type = function(index_type, adj, y, nperm = 100, max_dist=30, compartments=NULL, verbose=TRUE) {
    if (verbose) message(index_type)
    types = unique(y)
    i_index = which(y == index_type)
    i_shuffle = setdiff(seq_len(length(y)), i_index)
    X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+y) %>% as.matrix()
    colnames(X) = gsub('^y', '', colnames(X))
    freq = (colSums(X) / nrow(X))[types]
    freq_perm = map(seq_len(nperm), function(i) {
        set.seed(i)
        yperm = y
        if (is.null(compartments)) {
            yperm[i_shuffle] = sample(y[i_shuffle])
        } else {
            ## shuffle inside compartments, to preserve total composition within compartment
            .x = split(i_shuffle, compartments[i_shuffle]) %>%
                map(function(.i) {
                    ## CAUTION: if .i is a single number, sample will interpret it as 1:.i
                    if (length(.i) == 1) {
                        message('No shuffling is taking place, check code')
                        res = .i
                    } else {
                        res = sample(.i) ## shuffle non-index cells inside hub
                    }
                    names(res) = .i
                    return(res)
                }) %>%
                reduce(c)
            yperm[as.integer(names(.x))] <- y[.x]
        }
        X = adj[i_index, ] %*% Matrix::sparse.model.matrix(~0+yperm) %>% as.matrix() #%>% prop.table(1)
        colnames(X) = gsub('^yperm', '', colnames(X))
        (colSums(X) / nrow(X))[types]
    }) %>%
        purrr::reduce(rbind2)
    stats = tibble(
        type = types,
        freq,
        zscore = (freq - apply(freq_perm, 2, mean)) / apply(freq_perm, 2, sd),
        pval = exp(log(2) + (pnorm(-abs(zscore), log.p = TRUE, lower.tail = TRUE))), ## one-tailed
        fdr = p.adjust(pval)
    ) %>%
        cbind(dplyr::rename(data.frame(t(apply(freq_perm, 2, quantile, c(.025, .975)))), q025 = `X2.5.`, q975 = `X97.5.`)) %>% ## 95% CI
        subset(type != index_type) %>%
        dplyr::mutate(index_type = index_type) %>%
        dplyr::select(index_type, type, everything()) %>%
        arrange(fdr)
    return(stats)
}


coloc_all_types <- function(index_types, coords, y, nperm = 100, nsteps=1, max_dist=30, compartments=NULL, parallel=TRUE, verbose=TRUE) {
    if (parallel & length(index_types) > 1) {
        plan(multicore)
    } else {
        plan(sequential)
    }
    ## Define neighbors
    ## NOTE: max_dist only refers to directly adjacent neighbors
    adj = spatula::getSpatialNeighbors(coords, return_weights = TRUE)
    adj@x[adj@x > max_dist] = 0
    adj = Matrix::drop0(adj)
    adj@x = rep(1, length(adj@x))
    ## If nsteps>1, consider not only your adjacent neighbors
    ##   but also your neighbor’s neighbors etc.
    if (nsteps > 1) {
        adj = adj + Matrix::Diagonal(n = nrow(adj)) ## add self
        for (iter in seq_len(nsteps - 1)) {
            adj = adj %*% adj
        }
        ## Ignore weights. Only care if cell is a neighbor or not
        adj@x = rep(1, length(adj@x))
        ## Remove self as neighbor
        adj = adj - Matrix::Diagonal(n = nrow(adj))
        adj = Matrix::drop0(adj)
    }
    index_types %>%
        future_map(coloc_one_type, adj, y, nperm, max_dist, compartments, verbose, .options = furrr::furrr_options(seed = 1)) %>%
        rbindlist()  %>%
        identity
}



#analysis script

library(devtools)
library(spatula)
library(furrr)
 
all_x_y <- rbindlist(data_x_y)
 
#need to add a column to all_x_y with the ident for each cell. in this case mine is called 'named'
 
coloc_res_coarse<- coloc_all_types(
        index_type = unique(all_x_y$named),
        coords = all_x_y[, c("Centroid.X.µm", "Centroid.Y.µm")],
        y = all_x_y$named,
        compartments = NULL,
        max_dist = 40,
        nperm = 1000,
        parallel = TRUE
    )
 
library(data.table)
library(splitstackshape)
library(ggpubr)
library(rstatix)
library(DescTools)
 
plt_df<-coloc_res_coarse %>%
    subset(pval < 0.05) %>%
    dplyr::select(index_type, type, zscore) %>%
    spread(type, zscore, fill = 0) %>%
    column_to_rownames('index_type') %>%
    as.matrix
 
plt_df %>%
    Heatmap(col = circlize::colorRamp2(c(-20, 0, 20), c('blue', 'white', 'red')), border = T, cluster_columns = T, cluster_rows = T)
 
