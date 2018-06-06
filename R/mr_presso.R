mr_presso <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, OUTLIERtest = FALSE, DISTORTIONtest = FALSE, SignifThreshold = 0.05, NbDistribution = 1000, seed = NULL){

	if(!is.null(seed))
		set.seed(seed)

	if(SignifThreshold > 1)
		stop("The significance threshold cannot be greater than 1")

	if(length(BetaExposure) != length(SdExposure))
		stop("BetaExposure and SdExposure must have the same number of elements")

	# Functions
	"%^%" <- function(x, n) with(eigen(x), vectors %*% (values^n * t(vectors)))

	getRSS_LOO <- function(BetaOutcome, BetaExposure, data, returnIV){
		dataW <- data[, c(BetaOutcome, BetaExposure)] * sqrt(data[, "Weights"])
		X <- as.matrix(dataW[ , BetaExposure])
		Y <- as.matrix(dataW[ , BetaOutcome])
		CausalEstimate_LOO <- sapply(1:nrow(dataW), function(i) {
			(t(X[-i, ]) %*% X[-i, ])%^%(-1) %*% t(X[-i, ]) %*% Y[-i, ]
		})

		if(length(BetaExposure) == 1)
			RSS <- sum((Y - CausalEstimate_LOO * X)^2, na.rm = TRUE)
		else
			RSS <- sum((Y - rowSums(t(CausalEstimate_LOO) * X))^2, na.rm = TRUE)

		if(returnIV)
			RSS <- list(RSS, CausalEstimate_LOO)
		return(RSS)
	}


	getRandomData <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data){
		mod_IVW <- lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)), weights = Weights, data = data)
		dataRandom <- cbind(eval(parse(text = paste0("cbind(", paste0("rnorm(nrow(data), data[, \'", BetaExposure, "\'], data[ ,\'", SdExposure, "\'])", collapse = ", "), ", rnorm(nrow(data), fitted(mod_IVW), data[ ,\'", SdOutcome,"\']))"))), data$Weights)
		colnames(dataRandom) <- c(BetaExposure, BetaOutcome, "Weights")
		return(dataRandom)
	}

	# 0- Transforming the data + checking number of observations
	data <- data[, c(BetaOutcome, BetaExposure, SdOutcome, SdExposure)]
	data <- data[rowSums(is.na(data)) == 0, ]
	data[, c(BetaOutcome, BetaExposure)] <- data[, c(BetaOutcome, BetaExposure)] * sign(data[, BetaExposure[1]])
	data$Weights <- data$Weights <- 1/data[, SdOutcome]^2

	if(nrow(data) <= length(BetaExposure) + 2)
		stop("Not enough intrumental variables")

	if(nrow(data) >= NbDistribution)
		stop("Not enough elements to compute empirical P-values, increase NbDistribution")

	# 1- Computing the observed residual sum of squares (RSS)
	RSSobs <- getRSS_LOO(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, data = data, returnIV = OUTLIERtest)

	# 2- Computing the distribtion of expected residual sum of squares (RSS)
	randomData <- replicate(NbDistribution, getRandomData(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, SdOutcome = SdOutcome, SdExposure = SdExposure, data = data), simplify = FALSE)
	RSSexp <- sapply(randomData, getRSS_LOO, BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, returnIV = OUTLIERtest)
	if(OUTLIERtest)
		GlobalTest <- list(RSSobs = RSSobs[[1]], Pvalue = sum(RSSexp[1, ] > RSSobs[[1]])/NbDistribution)
	else
		GlobalTest <- list(RSSobs = RSSobs[[1]], Pvalue = sum(RSSexp > RSSobs[[1]])/NbDistribution)

	# 3- Computing the single IV outlier test
	if(GlobalTest$Pvalue < SignifThreshold & OUTLIERtest){
		OutlierTest <- do.call("rbind", lapply(1:nrow(data), function(SNV){
			randomSNP <- do.call("rbind", lapply(randomData, function(mat) mat[SNV, ]))
			if(length(BetaExposure) == 1){
				Dif <- data[SNV, BetaOutcome] - data[SNV, BetaExposure] * RSSobs[[2]][SNV]
				Exp <- randomSNP[, BetaOutcome] - randomSNP[, BetaExposure] * RSSobs[[2]][SNV]
			} else {
				Dif <- data[SNV, BetaOutcome] - sum(data[SNV, BetaExposure] * RSSobs[[2]][, SNV])
				Exp <- randomSNP[, BetaOutcome] - rowSums(randomSNP[, BetaExposure] * RSSobs[[2]][, SNV])
			}
			pval <- sum(Exp^2 > Dif^2)/length(randomData)
			pval <- cbind.data.frame(RSSobs = Dif^2, Pvalue = pval)
			return(pval)
		}))
		row.names(OutlierTest) <- row.names(data)
		OutlierTest$Pvalue <- apply(cbind(OutlierTest$Pvalue*nrow(data), 1), 1, min) # Bonferroni correction
	} else{
		OUTLIERtest <- FALSE
	}

	# 4- Computing the test of the distortion of the causal estimate
	mod_all <- lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data)
	if(DISTORTIONtest & OUTLIERtest){
		getRandomBias <- function(BetaOutcome, BetaExposure, SdOutcome, SdExposure, data, refOutlier){
			indices <- c(refOutlier, replicate(nrow(data)-length(refOutlier), sample(setdiff(1:nrow(data), refOutlier))[1]))
			mod_random <- lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", paste(BetaExposure, collapse = "+"))), weights = Weights, data = data[indices[1:(length(indices) - length(refOutlier))], ])
			return(mod_random$coefficients[BetaExposure])
		}
		refOutlier <- which(OutlierTest$Pvalue <= SignifThreshold)

		if(length(refOutlier) > 0){
			BiasExp <- replicate(NbDistribution, getRandomBias(BetaOutcome = BetaOutcome, BetaExposure = BetaExposure, data = data, refOutlier = refOutlier), simplify = FALSE)
			BiasExp <- do.call("rbind", BiasExp)

			mod_noOutliers <- lm(as.formula(paste0(BetaOutcome, " ~ -1 + ", BetaExposure)), weights = Weights, data = data[-refOutlier, ])
			BiasObs <- (mod_all$coefficients[BetaExposure] - mod_noOutliers$coefficients[BetaExposure]) / abs(mod_noOutliers$coefficients[BetaExposure])
			BiasExp <- (mod_all$coefficients[BetaExposure] - BiasExp) / abs(BiasExp)
			BiasTest <- list(`Outliers Indices` = refOutlier, `Distortion Coefficient` = 100*BiasObs, Pvalue = sum(abs(BiasExp) > abs(BiasObs))/NbDistribution)
		} else{
			BiasTest <- list(`Outliers Indices` = "No significant outliers", `Distortion Coefficient` = NA, Pvalue = NA)
		}
	}

	# 5- Formatting the results
	GlobalTest$Pvalue <- ifelse(GlobalTest$Pvalue == 0, paste0("<", 1/NbDistribution), GlobalTest$Pvalue)
	if(OUTLIERtest){
		OutlierTest$Pvalue <- replace(OutlierTest$Pvalue, OutlierTest$Pvalue == 0, paste0("<", nrow(data)/NbDistribution))
		if(DISTORTIONtest){
			BiasTest$Pvalue <- ifelse(BiasTest$Pvalue == 0, paste0("<", 1/NbDistribution), BiasTest$Pvalue)
			res <- list(`Global Test` = GlobalTest, `Outlier Test` = OutlierTest, `Distortion Test` = BiasTest)
		} else {
			res <- list(`Global Test` = GlobalTest, `Outlier Test` = OutlierTest)
		}
		if(nrow(data)/NbDistribution > SignifThreshold)
		warning(paste0("Outlier test unstable. The significance threshold of ", SignifThreshold, " for the outlier test is not achievable with only ", NbDistribution, " to compute the null distribution. The current precision is <", nrow(data)/NbDistribution, ". Increase NbDistribution."))
	} else {
		res <- GlobalTest
	}

	OriginalMR <- cbind.data.frame(BetaExposure, "Raw", summary(mod_all)$coefficients)
	colnames(OriginalMR) <- c("Exposure", "MR Analysis", "Causal Estimate", "Sd", "T-stat", "P-value")
	if(exists("mod_noOutliers"))
		OutlierCorrectedMR <- cbind.data.frame(BetaExposure, "Outlier-corrected", summary(mod_noOutliers)$coefficients)
	else
		OutlierCorrectedMR <- cbind.data.frame(BetaExposure, "Outlier-corrected", t(rep(NA, 4)))
	colnames(OutlierCorrectedMR) <- colnames(OriginalMR)
	MR <- rbind.data.frame(OriginalMR, OutlierCorrectedMR)
	row.names(MR) <- NULL

	res <- list(`Main MR results` = MR, `MR-PRESSO results` = res)
	return(res)
}
