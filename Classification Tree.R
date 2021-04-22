library(data.tree)
library(readxl)

# Function for plotting given tree
plot_tree<-function(TREE, N = TRUE, enumerate = FALSE, rounding = NULL){
  # TREE - tree to be plotted
  # N - if true, the function displays the number of observations in each node
  # enumerate - if true, the function displays numbers on leaves
  # rounding - number of decimal places that dependent values are rounded to
  if(is.null(TREE$plotted)){
    env <- new.env()
    assign('i', 1, envir=env)
    FUN <-function(x){
      SetEdgeStyle(x, label = "YES", node = GetAttribute(x, "YES"))
      SetEdgeStyle(x, label = "NO", node = GetAttribute(x, "NO"))
      if(enumerate & x$isLeaf){
        i <- get('i', envir=env)
        paste_enum <- paste0("Leaf ",i,".","\n")
        assign('i', i+1, envir=env)
        x$leaf <- i
      }else{
        paste_enum <- ""
      }
      if(!is.null(x$variable)){
        paste_cond <- paste0(x$variable," <= ",x$constraint, "?")
      }else{
        paste_cond <- ""
      }
      if(N){
        paste_N <- paste0("N = ", x$N,"\n")
      }else{
        paste_N <- "\n"
      }
      if(is.null(rounding)){
        new_name <- paste0(paste_N, "Y = ", x$y, 
                           "\n", paste_cond, paste_enum) 
      }else{
        new_name <- paste0(paste_N, "Y = ", round(x$y, digits = rounding), 
                           "\n", paste_cond, paste_enum) 
      }

      x$name <- new_name
    }
    TREE$Do(FUN)
  }
  TREE$plotted <- T
  plot(TREE)
}

# function for selecting the best split threshold for given attribute
find_threshold <- function(DATA, COL, DEP, minimum = 0){
  # DATA - data frame with variables COL and DEP
  # COL - independent variable for which the best split should be chosen
  # DEP - dependent variable
  # minimum - minimal number of observations in each node after the split
  DATA_ordered <- DATA[order(GetAttribute(DATA, COL)),]
  threshold <- GetAttribute(DATA_ordered, COL)[1] - 1
  mean_DEP <- mean(GetAttribute(DATA_ordered, DEP))
  RSS <- sum((GetAttribute(DATA_ordered, DEP)-mean_DEP)^2)
  RIGHT_N <- 0
  LEFT_N <- 0
  LEFT_mean <- 0
  RIGHT_mean <- 0
  if(2*minimum<nrow(DATA_ordered)){
    for(x in unique(GetAttribute(DATA_ordered, COL)[minimum:(nrow(DATA_ordered)-minimum)])){
      left <- DATA_ordered[GetAttribute(DATA_ordered, COL) <= x,]
      right <- DATA_ordered[GetAttribute(DATA_ordered, COL) > x,]
      left_mean <- mean(GetAttribute(left, DEP))
      left_sum <- sum((GetAttribute(left, DEP)-left_mean)^2)
      right_mean <- mean(GetAttribute(right, DEP))
      right_sum <- sum((GetAttribute(right, DEP)-right_mean)^2)
      RSS_partition <- left_sum+right_sum
      if((RSS_partition < RSS) & (nrow(left) >= minimum) & (nrow(right)) >= minimum){
        RSS <- RSS_partition
        threshold <- x
        LEFT_N <- nrow(left)
        RIGHT_N <- nrow(right)
        LEFT_mean <- mean(GetAttribute(left, DEP))
        RIGHT_mean <- mean(GetAttribute(right, DEP))
      }
    } 
  }
  ret <- list(threshold, LEFT_N, RIGHT_N, LEFT_mean, RIGHT_mean, RSS)
  names(ret) <- c('threshold', 'LEFT_N', 'RIGHT_N', 'LEFT_mean', 'RIGHT_mean', 'RSS')
  return(ret)
}

# function for selecting the best attribute and treshdold for partition
select_partition <- function(DATA, DEP, COLS = NULL, minimum = 0, alpha = 0){
  # DATA - data frame with independent variable and dependents variables
  # DEP - dependent variable
  # COLS  - vector of characters with attribute names (if want to use only a subset of predictors)
  # minimum - minimum number of observations in each leaf (if minimum < 1, minimal proportion of observations in the leaf)
  # alpha - weight in the cost function
  if((0 < minimum) & (minimum < 1)){
    minimum <- round(minimum*nrow(DATA))
  }
  if(is.null(COLS)){
    COLS <- colnames(DATA)
  }
  mean_DEP <- mean(GetAttribute(DATA, DEP))
  RSS <- sum((GetAttribute(DATA, DEP)-mean_DEP)^2)
  ret <- NULL
  for(x in COLS){
    if(x != DEP){
      summ <- find_threshold(DATA = DATA, COL = x, DEP = DEP, minimum = minimum)
      if((summ$LEFT_N > 0) & ((summ$RSS + alpha) < RSS)){
        RSS <- summ$RSS
        ret <- summ
        ret[7] <- x
      }
    }
  }
  if(!is.null(ret)){
    names(ret) <- c(names(ret)[1:6],"variable") 
  }
  return(ret)
}

# function for recursive sample partitioning and creating nodes
add_nodes <- function(TREE, DATA, DEP, COLS = NULL, minimum = 0, alpha = 0, 
                      depth = 1, max_depth = 100){
  # TREE - current node
  # DATA - data frame with variable DEP
  # DEP - predicted variable
  # COL - optional argument for providng a subset of predictors
  # minimum - minimal number of observations in each leaf or proportion of sample size in each leaf
  # alpha - weight in the cost function
  # depth - level of the current node
  # max_depth - maximal depth of the tree
    X <- select_partition(DATA = DATA, DEP = DEP, COLS = COLS, minimum = minimum, alpha = alpha)
    if(!is.null(X) & (depth < max_depth)){
      Left <- DATA[GetAttribute(DATA, X$variable) <= X$threshold,]
      Right <- DATA[GetAttribute(DATA, X$variable) > X$threshold,]
      TREE$variable <- X$variable
      TREE$constraint <- X$threshold
      TREE_LEFT <- TREE$AddChild("YES", N = nrow(Left), y = X$LEFT_mean)
      TREE_RIGHT <- TREE$AddChild("NO", N = nrow(Right), y = X$RIGHT_mean)
      add_nodes(TREE = TREE_LEFT, DATA = Left, DEP = DEP, COLS = COLS, 
                minimum = minimum, alpha = alpha, depth = (depth + 1), max_depth = max_depth)
      add_nodes(TREE = TREE_RIGHT, DATA = Right, DEP = DEP, COLS = COLS, 
                minimum = minimum, alpha = alpha, depth = (depth + 1), max_depth = max_depth)
    }
  return(TREE)
}

# function that creates the root and call recursive function for creating the whole tree
build_tree <- function(DATA, DEP, COLS = NULL, minimum = 0, alpha = 0, max_depth = 100){
  # data frame with dependent variable DEP
  # DEP - dependent variable
  # COL - optional argument for selecting a subset of predictors
  # minimum - minimal number of observations in each leaf or proportion of sample size in each leaf
  # alpha - weight for the cost function
  # max_depth - maximal depth of the tree
  if((0<minimum)&(minimum<1)){
    minimum <- round(nrow(DATA)*minimum, digits = 0)
  }
  TREE <- Node$new("Root", N = nrow(DATA), y = mean(GetAttribute(DATA, DEP)))
  add_nodes(TREE = TREE, DATA = DATA, DEP = DEP, COLS = COLS, 
            minimum = minimum, alpha = alpha, depth = 1, max_depth = max_depth)
  return(TREE)
}

# function for predicting the value of the dependent variable for given observation
predict <- function(TREE, obs){
  # TREE - regression tree to predict from
  # obs - single observation with necessary attributes
  while(!isLeaf(TREE)){
    children <- TREE$children
    if(GetAttribute(obs, TREE$variable) <= TREE$constraint){
      TREE <- children[[1]]
    }else{
      TREE <- children[[2]]
    }
  }
  if(!is.null(TREE$leaf)){
    cat(paste0("Leaf ",TREE$leaf,".","\n"))
  }
  return(TREE$y)
}
