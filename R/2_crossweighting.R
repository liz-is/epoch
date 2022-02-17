#' Perform crossweighting
#'
#' @param grnTab GRN dataframe, the result of running reconstructGRN or reconstructGRN_GENIE3
#' @param expDat genes-by-cells expression matrix
#' @param xdyn result of running findDynGenes
#' @param lag lag window on which to run cross-correlation. Cross-correlaiton computed from -lag to +lag.
#' @param min minimum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) less than min will not be negatively weighted.
#' @param max maximum of weighting window. Edges with offsets (or absolute offsets if symmetric_filter=TRUE) greater than max will have weights set to 0.
#' @param symmetric_filter whether or not to employ a symmetric weight scheme. If true, absolute offset is used in place of offset.
#' @param filter_thresh after crossweighting, edges with weights less than filter_thresh will be set to 0.
#'
#' @return grnTab with offset and weighted_score added
#' 
#' @export
#' @importFrom furrr future_pmap_dbl future_map2_dbl
#' @importFrom progressr progressor
#'
crossweight<-function(grnTab,
					expDat,
					xdyn,
					lag=floor(ncol(expDat)/5),
					min=ceiling(ncol(expDat)/50),
					max=floor(ncol(expDat)/12),
					symmetric_filter=FALSE,
					filter_thresh=0){

	# order expDat
	expDat<-expDat[,rownames(xdyn$cells)]

	grnTab$TG<-as.character(grnTab$TG)
	grnTab$TF<-as.character(grnTab$TF)
	
	TF_data <- lapply(grnTab$TF, function(TF) expDat[TF,])
	TG_data <- lapply(grnTab$TG, function(TG) expDat[TG,])
	
	message("Calculating offsets...")
	p <- progressr::progressor(steps = length(TF_data))
	offset <- furrr::future_map2_dbl(TF_data, TG_data, function(TF, TG){
	  p()
	  cross_corr(TF, TG, lag = lag)
	  })
	grnTab$offset <- offset
	
	message("Calculating weights...")
	p <- progressr::progressor(steps = nrow(grnTab))
	weighted_score <- furrr::future_pmap_dbl(grnTab, function(zscore, offset, ...) {
	  p()
	  score_offset(score = zscore, offset = offset, min = min, max = max, 
	               symmetric_filter = symmetric_filter)
	  })
  grnTab$weighted_score <- weighted_score
  
	grnTab<-grnTab[grnTab$weighted_score>filter_thresh,]

	return(grnTab)
}


cross_corr<-function(TF, TG, expDat,lag, ...){
	x<-ccf(TF,TG,lag,pl=FALSE)

	df<-data.frame(lag=x$lag,cor=abs(x$acf))
  df<-df[order(df$cor,decreasing=TRUE),]
  offset<-mean(df$lag[1:ceiling((2/3)*lag)])

  return(offset)
}

score_offset<-function(score,offset,min=2,max=20,symmetric_filter=FALSE){

	if(symmetric_filter){
		offset<-abs(offset)
	}

	if (offset<=min){
		res<-score
	}else if (offset>=max){
		res<-0
	}else{
		# linear weighting scheme according to y=(-x/(max-min))+1
		weight<-(-offset/(max-min))+1
		res<-score*weight
	}

	res
}


# estimates min and max values for crossweighting
# for now assumes uniform cell density across pseudotime/only considers early time
# this needs to be refined if it's to be useful...
#' @export
crossweight_params<-function(expDat,xdyn,pseudotime_min=0.005,pseudotime_max=0.01){

	expDat<-expDat[,rownames(xdyn$cells)]
	ncells<-nrow(xdyn$cells)
	min<-nrow(xdyn$cells[xdyn$cells$pseudotime<pseudotime_min,])
	max<-nrow(xdyn$cells[xdyn$cells$pseudotime<pseudotime_max,])

	params<-list(min=min,max=max)

}







