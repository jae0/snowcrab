
  logbook.prorate = function(lgbk) {
    
    # manual pro-rate:
  
    # fix problems with slip weights
      new.slip.sums = tapply2( x=lgbk, indices=c("doc_id"), var="slip_weight_lbs", newvars=c("new.slip.wgt"),
                       func=function(x){ sum(unique(x)) } )

      lgbk = merge( x=lgbk, y=new.slip.sums, by="doc_id", sort=F, all.x=T, all.y=F)
      i = which(( lgbk$new.slip.wgt - lgbk$slip_weight_lbs) != 0 )
      j = which( duplicated( lgbk$log_efrt_std_info_id,  ))
      bad.data = intersect ( i, j )
      if (length(bad.data) > 0 ) lgbk = lgbk[-bad.data ,]

      lgbk$slip_weight_lbs = lgbk$new.slip.wgt
      lgbk$new.slip.wgt = NULL
  
      # compute rescaling factors that ensure sums are still correct after rescaling
      # counts
      pr.count = tapply2( x=lgbk, indices=c("doc_id"), var="log_efrt_std_info_id", newvars=c("pr.count"),
                       func=function(x){ length(unique(x)) } )
      # sums
      pr.sum = tapply2( x=lgbk, indices=c("doc_id"), var="est_weight_log_lbs", newvars=c("pr.sum"),
                       func=function(x){ sum(x, na.rm=F) } )

      # merge and rescale
      pr = merge(x=pr.count, y=pr.sum, by="doc_id", sort=F, all=F)
      pr$doc_id = as.character(pr$doc_id )

      lgbk = merge(x=lgbk, y=pr, by="doc_id", all.x=T, all.y=F, sort=F)
      lgbk$pr.fraction = lgbk$est_weight_log_lbs / lgbk$pr.sum
      lgbk$pro_rated_slip_wt_lbs = lgbk$pr.fraction * lgbk$slip_weight_lbs

      # second pro-rating of data with missing values/nulls, etc
      na.logs = which (!is.finite ( lgbk$pr.sum ) )
      if (length( na.logs) > 0 ) {
        lgbk$pro_rated_slip_wt_lbs[ na.logs ] = lgbk$slip_weight_lbs[na.logs ] / lgbk$pr.count [  na.logs]
      }

      bad.usage.codes = which ( lgbk$catch_usage_code!=10 )
      if (length(bad.usage.codes) > 0 ) lgbk = lgbk[-bad.usage.codes, ]

    return( lgbk )
  }



