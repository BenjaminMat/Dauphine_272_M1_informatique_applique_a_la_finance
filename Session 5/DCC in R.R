library(rmgarch)

data_path <- '/Users/renzhi/Desktop/Python2025_2026/Renzhi Material/Session 4/Data'

Ret_assets <- read.csv2(file.path(data_path, "ret_asset.csv"), sep = ',')
Ret_assets$Date <- as.Date(Ret_assets$Date, format = "%Y-%m-%d")
Ret_assets[, 2:ncol(Ret_assets)] <- lapply(Ret_assets[, 2:ncol(Ret_assets)], as.numeric)

Ret_mkt <- read.csv2(file.path(data_path, "ret_mkt.csv"), sep = ',')
Ret_mkt$Date <- as.Date(Ret_mkt$Date, format = "%Y-%m-%d")
Ret_mkt[, 2] <- as.numeric(Ret_mkt[, 2])

cal_rho <- function(dfa, dfm, window = 250) {
  p <- ncol(dfa)
  n <- nrow(dfa)
  
  df_rho <- matrix(NA_real_, nrow = n, ncol = p)  # same shape as dfa
  
  # Univariate margins: sGARCH(1,1), mean=0, normal
  uni_spec <- ugarchspec(
    variance.model   = list(model = "sGARCH", garchOrder = c(1, 1)),
    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
    distribution.model = "norm"  
  )
  mspec <- multispec(replicate(2, uni_spec))
  
  # DCC(1,1)
  spec <- dccspec(uspec = mspec, dccOrder = c(1, 1), distribution = "mvt")
  
  for (col in 2:p) {
    data_matrix <- cbind(dfa[, col], dfm[, 2])  # asset vs market
    
    if (nrow(stats::na.omit(data_matrix)) <= window) next
    
    for (j in window:(nrow(data_matrix) - 1)) {
      idx <- (j - window + 1):j
      rolling_window <- data_matrix[idx, , drop = FALSE]
      if (any(!is.finite(rolling_window))) next
      
      target_row <- j + 1
      
      tryCatch({
        fit <- dccfit(spec, data = rolling_window)
        
        if (!inherits(fit, "DCCfit")) stop("Non-converged univariate fits (uGARCHmultifit).")
        
        fcast <- dccforecast(fit, n.ahead = 1, n.roll = 0)
        Hf    <- fcast@mforecast$H[[1]][,,1]  # 2x2 matrix: cov_t+1|t
        
        rho_forecast <- Hf[1, 2] / sqrt(Hf[1, 1] * Hf[2, 2])
        df_rho[target_row, col] <- rho_forecast
        
        message(sprintf("rho: row %d / col %d", target_row, col))
      }, error = function(e) {
        df_rho[target_row, col] <- NA_real_
        message(sprintf("Error at row %d col %d: %s", target_row, col, e$message))
      })
    }
  }
  
  df_rho
}

df_rho_mat <- cal_rho(Ret_assets, Ret_mkt, window = 250)

df_rho <- as.data.frame(df_rho_mat)
df_rho[, 1] <- Ret_assets[, 1]           # Date to first column
colnames(df_rho) <- colnames(Ret_assets)


write.csv(df_rho, file.path(data_path, "DCC_Rho.csv"), row.names = FALSE)




