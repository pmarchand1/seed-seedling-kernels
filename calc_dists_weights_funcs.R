# Functions to calculate trap-tree distances and 
#  Ripley edge correction weights

# For a circle of radius r centered on (0, 0),
#  get the anglular portion contained in the rectangle (0,0) to (x,y)
get_angle_in_rect <- function(r, x, y) {
    ifelse(r^2 >= x^2 + y^2, 0, asin(y / pmax(y, r)) - acos(x / pmax(x, r)))
}

# Calculate distance matrix and weights (fraction of radius in plot) 
#  for each trap-tree pair
calc_dists_weights <- function(traps, trees) {
    d2min <- 0.01
    
    dist_sq <- outer(traps$x, trees$x, "-")^2 + outer(traps$y, trees$y, "-")^2
    dist_sq[dist_sq < d2min] <- d2min
    
    r <- sqrt(dist_sq)
    
    xr <- xmax - traps$x
    xl <- traps$x - xmin
    yt <- ymax - traps$y
    yb <- traps$y - ymin
    rmax <- sqrt(pmax(xr^2 + yt^2, xl^2 + yt^2, xr^2 + yb^2, xl^2 + yb^2))
    
    ntree <- nrow(trees)
    xr <- matrix(rep(xr, ntree), ncol = ntree)
    xl <- matrix(rep(xl, ntree), ncol = ntree) 
    yt <- matrix(rep(yt, ntree), ncol = ntree)
    yb <- matrix(rep(yb, ntree), ncol = ntree)
    
    inv_wgt <- 0.5/pi * (get_angle_in_rect(r, xr, yt) + get_angle_in_rect(r, xl, yt)
                         + get_angle_in_rect(r, xr, yb) + get_angle_in_rect(r, xl, yb))
    list(r = r, rmax = rmax, wgt = 1 / inv_wgt)
}