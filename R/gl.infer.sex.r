#'@name gl.infer.sex
#'@title Uses sex-linked loci to infer sex of individuals
#'@description
#' This function uses the output of function gl.keep.sexlinked (list of 5 
#' elements) to infer the sex of all individuals in the dataset. 
#' It uses 3 types of sex-linked loci (W-/Y-linked, Z-/X-linked, and 
#' gametologs), assigns a preliminary genetic sex for each type of sex-linked
#' loci available, and outputs an agreed sex.
#' 
#' This function produces as output a dataframe with individuals in rows and 11 
#' columns.
#'
#' @param gl_sexlinked The output of function gl.keep.sexlinked (complete 
#' list with 5 elements). See explanation in "Details" section [required].
#' @param system String that declares the sex-determination system of the 
#' species: 'zw' or 'xy' [required].
#' @param seed User-defined integer for repeatability purposes. If not provided
#' by user, it is chosen randomly by default. See "Details" section.
#'
#' @details
#' Parameter \code{gl_sexlinked} must be the name of the output object (a 
#' list of 5 elements) produced by function \code{gl.keep.sexlinked}. Parameter 
#' \code{seed} must be an integer that will be used on the KMeans algorithm used
#' by the function. We highly recommend choosing the seed to guarantee 
#' repeatability.
#' 
#' Note that this function was created with the explicit intent that a human 
#' checks the evidence for the sex assignments that do NOT agree for all 
#' types of sex-linked loci (called "indefinite sex assignments" and denoted 
#' as "*M" or "*F" in the last column of dataframe output). This human can then 
#' use their criterion to validate these assignments.
#' 
#'
#'\strong{ Function's output }
#'
#' This function creates a dataframe with one row per individual and 11 columns: \itemize{
#' \item {id > Individuals' ID.}
#' 
#' \item {w.linked.sex or y.linked.sex > Sex inferred using w-linked or y-linked
#' loci.}
#' \item {#called > Number of W-linked or Y-linked loci for which the individual 
#' had a called genotype (cf. missing genotype).}
#' \item {#missing > Number of W-linked or Y-linked loci for which the 
#' individual had a missing genotype (cf. called genotype).} 
#' 
#' \item {z.linked.sex or x.linked.sex > Sex inferred using z-linked or x-linked
#' loci.}
#' \item {#Hom.z or #Hom.x > Number of z-linked or x-linked loci for which the 
#' individual is homozygous.}
#' \item {#Het.z or #Het.x > Number of z-linked or x-linked loci for which the 
#' individual is heterozygous.}
#' 
#' \item {gametolog.sex > Sex inferred using gametologs.} 
#' \item {#Hom.g > Number of gametologous loci for which the individual is 
#' homozygous.}
#' \item {#Het.g > Number of gametologous loci for which the individual is 
#' heterozygous.}
#' 
#' \item {agreed.sex > Agreed sex: 'F' or 'M' if all preliminary sex-assignments
#' match (i.e., definite sex assignment), and '*F' or '*M' if NOT all 
#' preliminary sex-assignments match (i.e., indefinite sex assignment).} 
#' }
#' 
#' @return A dataframe.
#' @author Custodian: Diana Robledo-Ruiz  -- Post to
#'   \url{https://groups.google.com/d/forum/dartr}
#'
#' @examples
#' LBP_sexLinked <- gl.keep.sexlinked(x = LBP, system = "xy", plot.display = TRUE, ncores = 1)
#' inferred.sexes <- gl.infer.sex(gl_sexlinked = LBP_sexLinked, system = "xy", seed = 100)
#' inferred.sexes
#' 
#' @references
#' \itemize{
#' \item Robledo‐Ruiz, D. A., Austin, L., Amos, J. N., Castrejón‐Figueroa, J.,
#'  Harley, D. K., Magrath, M. J., Sunnucks, P., & Pavlova, A. (2023). 
#'  Easy‐to‐use R functions to separate reduced‐representation genomic datasets
#'   into sex‐linked and autosomal loci, and conduct sex assignment. Molecular 
#'   Ecology Resources, 00, 1-21.
#'  }
#' 
#' @importFrom stats kmeans
#' @importFrom stats na.omit
#' 
#' @export
gl.infer.sex <- function(gl_sexlinked, 
                         system = NULL, 
                         seed = NULL) {
  
  # Parameters check
  if(is.null(system)){
    stop("You must specify the sex-determination system with the parameter 'system' ('zw' or 'xy').")
  } else {
    if(!(system == 'zw' | system == 'xy')){
      stop("Parameter 'system' must be 'zw' or 'xy'.")
    }
  }
  
  # Random seed if not specified by user
  if(is.null(seed)) {
    seed <- sample.int(65535, 1)
  }
  
  if(system == "xy") {
    gl1 <- gl_sexlinked$y.linked
    gl2 <- gl_sexlinked$x.linked
  }
  
  if(system == "zw") {
    gl1 <- gl_sexlinked$w.linked
    gl2 <- gl_sexlinked$z.linked
  }
  
  # Gametologs
  gl3    <- gl_sexlinked$gametolog
  table  <- gl_sexlinked$results.table  # Retrieve table
  
  
  if(system == "zw") {
    all    <- table[table$zw.gametolog == TRUE, ]  # Gametologs
  } else {
    all    <- table[table$xy.gametolog == TRUE, ]  # Gametologs
  }
  all    <- all[order(all$stat.p.adjusted), ]      # Order from smallest p-value
  useful <- row.names(all[1:5, ])                  # Keep name of only top 5 gametologs
  
  
  # Make sex assignment per type of sex-linked loci (Functions declared below)
  # W/Y-linked
  if (gl1@n.loc >= 1){
    w <- W.sex(gl1, system = system)
  } else {
    message("Not enough W-linked/Y-linked loci (need at least 1). Assigning NA...")
    w <- data.frame(W.sex = rep(NA, length(gl1@ind.names)),
                    n0.w  = rep(NA, length(gl1@ind.names)),
                    n1.w  = rep(NA, length(gl1@ind.names)))
  }
  # Z/X-linked
  if (gl2@n.loc >= 2){
    z <- Z.sex(gl2, system = system, seed = seed)
  } else {
    message("Not enough Z-linked/X-linked loci (need at least 2). Assigning NA...")
    z <- data.frame(Z.sex = rep(NA, length(gl2@ind.names)),
                    n1.z  = rep(NA, length(gl2@ind.names)),
                    n0.z  = rep(NA, length(gl2@ind.names)))
  }
  # Gametologs
  if (gl3@n.loc >= 5){
    g <- g.sex(gl3, system = system, seed = seed, useful = useful)
  } else {
    message("Not enough gametologs (need at least 5). Assigning NA...")
    g <- data.frame(g.sex = rep(NA, length(gl3@ind.names)),
                    n1.g  = rep(NA, length(gl3@ind.names)),
                    n0.g  = rep(NA, length(gl3@ind.names)))
  }
  
  # Put them all together
  A <- data.frame(w, z, g)
  
  # Function to conciliate assignments
  Fun <- function(x, y, z){
    d  <- c(x, y, z)
    yy <- data.frame(x = c("F", "M"),
                     freq = c(sum(d == "F", na.rm = TRUE),
                              sum(d == "M", na.rm = TRUE)))
    yy <- yy[order(-yy$freq), ]
    value <- if(length(unique(na.omit(d))) == 1) unique(na.omit(d)) else sprintf('*%s', yy[1, 1])
    return(value)
  }
  
  # Add last column of agreed sexes
  A$agreed.sex <- mapply(Fun, A$W.sex, A$Z.sex, A$g.sex)
  
  if(system == 'xy'){
    names <- c('y.linked.sex',  '#called', '#missing',
               'x.linked.sex',  '#Het.x',   '#Hom.x',
               'gametolog.sex', '#Het.g',   '#Hom.g', 'agreed.sex')
  } else {
    names <- c('w.linked.sex',  '#called', '#missing',
               'z.linked.sex',  '#Het.z',  '#Hom.z',
               'gametolog.sex', '#Het.g',  '#Hom.g', 'agreed.sex')
  }
  
  colnames(A) <- names
  
  A <- cbind(row.names(A), A)
  
  colnames(A)[1] <- "id"
  
  message("***FINISHED***")
  
  return(A)
}





############################### 1. W.sex function
### Map NAs (missing) and scored (called) to 1Dim in [-1,1], IF x<0, F else M

W.sex <- function(gl, system = NULL){
  w <- as.matrix(gl)
  w[is.na(w)] <- 3
  
  n0.w <- rowSums(w == 0 | w == 2 | w == 1, na.rm = TRUE)
  n1.w <- rowSums(w == 3, na.rm = TRUE)
  
  # Calculate proportion
  sex.score <- function(f, m){
    return( (-f+m)/(f+m) )
  }
  
  c2 <- sex.score(n0.w, n1.w)
  
  if(system == 'xy'){
    lab0 <- 'M'
    lab1 <-'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }
  
  W.sex <- ifelse(c2 < 0, lab0, lab1)
  
  Y <- data.frame(W.sex, n0.w, n1.w)
  return(Y)
}

############################### 2. Z.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Hom

Z.sex <- function(gl, system = NULL, seed = 42){
  z <- as.matrix(gl)
  
  n0.z = rowSums(z == 0 | z == 2, na.rm = TRUE)
  n1.z = rowSums(z == 1, na.rm = TRUE)
  
  Z_unclean <- t(apply(data.frame(n0.z, n1.z), 1, function(x) x / sum(x) ))
  Z <- na.omit(Z_unclean)
  Zna <- Z_unclean[rowSums(is.na(Z_unclean)) > 0,]
  
  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z, 2)
  
  i <- names(which.max(n1.z)) # Largest proportion of '1'
  label <- km$cluster[i]
  
  if(system == 'xy'){
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }
  
  # Assign sex
  Z.sex <- ifelse( km$cluster  == label, lab1, lab0)
  
  # Assign NA to inds with NA
  if (!is.null(nrow(Zna)) && nrow(Zna) > 0 ){
    for (i in 1:nrow(Zna)) {
      Z.sex[length(Z.sex)+1] <- NA
      names(Z.sex)[length(names(Z.sex))] <- rownames(Zna)[i]
    }
  }
  
  # Fuse them
  Y <- data.frame(n1.z, n0.z)
  
  # Add NAs in appropriate row
  Y$Z.sex <- "STOP"
  for (i in 1:nrow(Y)){
    Y[i, "Z.sex"] <- Z.sex[rownames(Y)[i]]
  }
  
  Y <- Y[, c("Z.sex", "n1.z", "n0.z")]
  
  return(Y)
}

############################### 3. ZWg.sex function
### Map Hom and Het to 2Dim and apply kmeans. Choose the label from maximum Het

g.sex <-  function(gl, system = NULL, seed = 42, useful = useful) {
  
  z <- as.matrix(gl)
  z <- z[, useful]
  
  n0.g = rowSums(z == 0 | z == 2, na.rm = TRUE)
  n1.g = rowSums(z == 1, na.rm = TRUE)
  
  Z_unclean <- t(apply(data.frame(n0.g, n1.g), 1, function(x) x / sum(x) ))
  Z <- na.omit(Z_unclean)
  Zna <- Z_unclean[rowSums(is.na(Z_unclean)) > 0,]
  
  # Apply k-means
  set.seed(seed)
  km <- kmeans(Z, 2)
  
  i <- names(which.max(n1.g)) # Largest proportion of '1'
  label <- km$cluster[i]
  
  if(system == 'xy'){
    lab0 <- 'M'
    lab1 <- 'F'
  } else {
    lab0 <- 'F'
    lab1 <- 'M'
  }
  
  # Assign sex (HERE IS THE OPPOSITE)
  Z.sex <- ifelse( km$cluster  == label, lab0, lab1)
  
  # Assign NA to inds with NA
  if (!is.null(nrow(Zna)) && nrow(Zna) > 0 ){
    for (i in 1:nrow(Zna)) {
      Z.sex[length(Z.sex)+1] <- NA
      names(Z.sex)[length(names(Z.sex))] <- rownames(Zna)[i]
    }
  }
  
  # Fuse them
  Y <- data.frame(n1.g, n0.g)
  
  # Add NAs in appropriate row
  Y$g.sex <- "STOP"
  for (i in 1:nrow(Y)){
    Y[i, "g.sex"] <- Z.sex[rownames(Y)[i]]
  }
  
  Y <- Y[, c("g.sex", "n1.g", "n0.g")]
  
  return(Y)
}