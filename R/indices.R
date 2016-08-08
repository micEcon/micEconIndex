micEconIndex <- function( prices, quantities, base, data, method, na.rm,
   weights, EKS = FALSE ) {

   if( length( prices ) != length( quantities ) ) {
      stop( "arguments 'prices' and 'quantities' must have the same length" )
   }

   checkNames( c( prices, quantities ), names( data ) )

   n <- length( prices )

   numerator <- numeric( nrow( data ) )
   denominator <- numeric( nrow( data ) )

   if( method %in% c( "Laspeyres", "Paasche" ) ) {
      for( i in 1:n ) {
         pt <- data[[ prices[ i ] ]]
         p0 <- mean( data[[ prices[ i ] ]][ base ], na.rm = na.rm )
         if( is.na( p0 ) ){
            warning( "base of '", prices[ i ], "' is NA" )
         }
         if( is.infinite( p0 ) ){
            warning( "base of '", prices[ i ], "' is infinite" )
         }
         qt <- data[[ quantities[ i ] ]]
         q0 <- mean( data[[ quantities[ i ] ]][ base ], na.rm = na.rm )
         if( is.na( q0 ) ){
            warning( "base of '", quantities[ i ], "' is NA" )
         }
         if( is.infinite( q0 ) ){
            warning( "base of '", quantities[ i ], "' is infinite" )
         }
         if( method == "Laspeyres" ) {
            if( is.na( q0 ) ) {
               numerator[ pt != 0 | is.na( pt ) ] <- NA
               if( p0 != 0 | is.na( p0 ) ) {
                  denominator <- NA
               }
            } else if( q0 != 0 ) {
               numerator <- numerator +  pt * q0
               denominator <- denominator + p0 * q0
            }
         } else if( method == "Paasche" ) {
            if( is.na( p0 ) ) {
               denominator[ qt != 0 | is.na( qt ) ] <- NA
            } else if( p0 != 0 ) {
               denominator <- denominator + p0 * qt
            }
            numerator[ is.na( pt ) & ( qt != 0 | is.na( qt ) ) ] <- NA
            selection <- ( pt != 0 & !is.na( pt ) ) & ( qt != 0 | is.na( qt ) )
            numerator[ selection ] <- numerator[ selection ] +
                  pt[ selection ] * qt[ selection ]
         }
      }
      result <- numerator / denominator
      if( weights ) {
         numerator[ is.na( denominator ) ] <- NA
         weightData <- data.frame( obsNo = c( 1:nrow( data ) ) )
         rownames( weightData ) <- rownames( data )
         for( i in 1:n ) {
            if( method == "Laspeyres" ) {
               weightData[[ prices[ i ] ]] <- data[[ prices[ i ] ]] *
                  mean( data[[ quantities[ i ] ]][ base ], na.rm = na.rm ) /
                  numerator
            } else if( method == "Paasche" ) {
               weightData[[ prices[ i ] ]] <- data[[ prices[ i ] ]] *
                  data[[ quantities[ i ] ]] / numerator
            }
         }
         weightData$obsNo <- NULL
         attributes( result )$weights <- weightData
      }
   } else if( method == "Fisher" & ! EKS ) {

      pL <- priceIndex( prices, quantities, base, data, method = "Laspeyres",
         na.rm = na.rm, weights )
      pP <- priceIndex( prices, quantities, base, data, method = "Paasche",
         na.rm = na.rm, weights )
      result <- sqrt( pL * pP )
      if( weights ) {
         attributes( result )$weights <-
            0.5 * attributes( pL )$weights +
            0.5 * attributes( pP )$weights
      }
    } else if( method == "Fisher" & EKS ) {
        result <- fisherEKS(quantities, prices, data)
        # Counterintuitively, quants and prices are reversed because
        # of the way quantityIndex() is called
   } else {
      stop( "argument 'method' must be either 'Laspeyres', 'Paasche'",
         " or 'Fisher'" )
   }
   names( result ) <- rownames( data )
   return( result )
}

priceIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE, weights = FALSE ) {

   checkNames( c( prices, quantities ), names( data ) )

   result <- micEconIndex( prices, quantities, base, data, method, na.rm,
      weights, EKS = FALSE   )

   return( result )
}

quantityIndex <- function( prices, quantities, base, data,
   method = "Laspeyres", na.rm = FALSE, weights = FALSE, EKS = FALSE ) {

   checkNames( c( prices, quantities ), names( data ) )

   result <- micEconIndex( quantities, prices, base, data, method, na.rm,
      weights, EKS = EKS )

   return( result )
}


# consolMatrix <- function(x) {
#   if (!is.data.table(x)) x <- as.data.table(x)
#   xRet <- x[, .(.N), by = names(x)]
#   nRet <- matrix(xRet$N, ncol = 1)
#   xRet[, N := NULL]
#   list(mat = as.matrix(xRet), freq = nRet)
# }



fisherEKS <- function( prices, quantities, data) {

    if (all(!duplicated(data[, prices])) & all(!duplicated(data[, quantities]))) {

      ret <- fisherEKSdense(
        as.matrix(data[, quantities]),
        as.matrix(data[, prices]),
        matrix(rep(1, nrow(data)), ncol = 1),
        matrix(rep(1, nrow(data)), ncol = 1),
        1:nrow(data) - 1, # zero indexing of C++
        1:nrow(data) - 1
      )

    } else {
      stop("TODO")
    }

   return(as.vector(ret))
}
