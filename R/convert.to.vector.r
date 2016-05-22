
    convert.to.vector = function(Q, nodes, type="default") {

      x = data.frame( id=as.vector(nodes))
      x$order = seq(1:length(nodes))

      x$id = as.character(x$id)

      if (type=="default")  {
        v = as.data.frame(as.table(Q))
        names(v) = c("r","c","n")
        v$id = paste(v$c, v$r, sep=".")
      } else if (type=="observer.data") {
        v = Q
        names(v) = c("yr", "id", "n")
      }

      x = merge(x, v[,c("id", "n")], by="id", all.x=T, all.y=F, sort=F)
      x = x[order(x$order),]
      rn = x$id
      x = as.vector(x$n)
      x[ !is.finite(x) ] = 0
      x[ x<0.0001 ] = 0
      names(x) = rn

      return(x)
    }


