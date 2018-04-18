HEALTHY = 0
INFECTEDA1 = 1
INFECTEDA2 = 2
DEAD = 3
EMPTY = 4
numberOfA2 = 5 # Initial Number of A2?
responseTime  = 4*7 #4 weeks
responseTimeGrid <- array(numeric(), c(0,0))

HIV = function(n, probHIV, probEmpty, probInfect, probReplace, t) {
  # FIRE simulation
  body  = initBody( n, probHIV, probEmpty)
  responseTimeGrid <<- body * responseTime
  
  grids = array(dim=c(n,n,t+1))

  
  grids[,,1] = body
  
  for (i in 2:(t+1)) {
    bodyExtended = periodicLat(body)
    responseTimeGrid <<- periodicLat(responseTimeGrid)
    body = applyExtended(bodyExtended, probInfect, probReplace, i-1)
    grids[,,i] = body
  }
  
  return(grids)
}

applyExtended = function(latExt, probInfect, probReplace, layer) {
  # APPLYEXTENDED - Function to apply 
  # spread() to every interior
  # site of square array latExt and to return the resulting array
  n = nrow(latExt) - 2
  newmat = matrix(c(rep(0,n*n)), nrow = n)
  for (j in 2:(n + 1)) {
    for (i in 2:(n + 1)) {
      site = latExt[i, j]
      N = latExt[i - 1, j]
      NE = latExt [i-1, j+1]
      NW = latExt[i-1, j-1]
      E = latExt[i, j + 1]
      S = latExt[i + 1, j]
      W = latExt[i, j - 1]
      SW = latExt[i+1, j-1]
      SE = latExt[i+1, j+1]
      newmat[i - 1, j - 1] = spread(site, N, NE, NW, E, S, W, SW, SE, probInfect, probReplace, i, j)
      if (layer>1){
        move(site, N, NE, NW, E, S, W, SW, SE, grids, i, j, layer)
      }
    }
  }
  return(newmat)
}

initBody = function( n, probHIV, probEmpty){
  #INITBODY returns an n-by-n grid of values?\
  # HEALTHY (healthy cell), 
  # INFECTEDA1 (infected cell), or INFECTEDA2 (infected cell), or DEAD (dead cell)
  # Pre:	n is the size (number of rows or columns) of the square grid and is positive.
  # 	probHIV is the probability that a cell is initially infected 
  # Post:	 A grid as described above was returned.
  
  cells = matrix(runif(n^2) < 1-probEmpty,nrow=n,ncol=n)
  infected = cells * matrix(runif(n^2) < probHIV,nrow=n,ncol=n)
  healthy = cells-infected
  empty = 1-cells
  # 1 where infected
  #infected = matrix(runif(n^2) < probHIV,nrow=n,ncol=n)
  # 1 where healthy
  #healthy = 1 - (infected + empty)

  # EMPTY, TREE, or BURNING
  body =  healthy * HEALTHY + infected * INFECTEDA1 + empty * EMPTY
  
  return (body)
}

periodicLat = function(lat) {
  # PERIODICLAT returns extended lattice
  extendRows = rbind(lat[nrow(lat),],lat,lat[1,])
  extlat = cbind(extendRows[,ncol(extendRows)],extendRows,extendRows[,1])
  return(extlat)
}

pointsForGrid = function(grid,val) {
  # Helper function for showGraphs
  # Pre: grid is a square grid
  #       val is a possible element value
  # Pre: The function has returnd a list of xcoords and ycoords,
  #      coordinates of cells that have value val.
  #      To match matrix values, the rows are reversed.
  xcoords = vector()
  ycoords = vector()
  for (row in 1:nrow(grid)) {
    for (col in 1:ncol(grid)) {
      if (grid[row,col] == val) {
        xcoords[length(xcoords)+1] = col
        ycoords[length(ycoords)+1] = nrow(grid) - row
      }
    }
  }
  return(list(xcoords,ycoords))
}

showGraphs = function(graphList, n) {
  # SHOWGRAPHS - Function to perform animation of grids in graphList
  m = dim(graphList)[3]
  for (k in 1:m) {
    g = graphList[,,k]
    healthy = pointsForGrid(g, HEALTHY)
    infecteda1 = pointsForGrid(g,INFECTEDA1)
    infecteda2 = pointsForGrid(g,INFECTEDA2)
    dead = pointsForGrid(g,DEAD)
    empty = pointsForGrid(g,EMPTY)
    plot(healthy[[1]],healthy[[2]],pch=19,col="green",
         xlim=c(0,n+1),ylim=c(0,n+1))
    points(infecteda1[[1]],infecteda1[[2]],col="red",pch=23,bg="orange")
    points(infecteda2[[1]],infecteda2[[2]],col="red",pch=23,bg="purple")
    points(dead[[1]],dead[[2]],col="red",pch=23,bg="black")
    points(empty[[1]],empty[[2]],col="black",pch=23,bg="white")
    Sys.sleep(0.2)
  }
}

move = function(site, N, NE, NW, E, S, W, SW, SE, grid, i, j, layer) {
  if(site == HEALTHY || site == INFECTEDA1 || site == INFECTEDA2 || site == DEAD) {
    space = c(0)
    # Finds which spaces around the site are EMPTY and stores them in a vector
    if(N == EMPTY) {space = c(space, N)}
    if(NE == EMPTY) {space = c(space, NE)}
    if(NW == EMPTY) {space = c(space, NW)}
    if(E == EMPTY) {space = c(space, E)}
    if(S == EMPTY) {space = c(space, S)}
    if(W == EMPTY) {space = c(space, W)}
    if(SE == EMPTY) {space = c(space, SE)}
    if(SW == EMPTY) {space = c(space, SW)}
    disp = space[sample(1:length(space),1)]
    print(i)
    print(j)
    if (disp == N) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i-1 , j]
      responseTimeGrid[i-1 , j] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i-1 , j,layer]
      grid[i-1 , j,layer] <- hold
    }
    if (disp == S) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i+1 , j]
      responseTimeGrid[i+1 , j] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i+1 , j,layer]
      grid[i+1 , j,layer] <- hold
    }
    if (disp == E) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i , j+1]
      responseTimeGrid[i , j+1] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i , j+1,layer]
      grid[i , j+1,layer] <- hold
    }
    if (disp == W) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i , j-1]
      responseTimeGrid[i , j-1] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i , j-1,layer]
      grid[i , j-1,layer] <- hold
    }
    if (disp == NE) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i-1 , j+1]
      responseTimeGrid[i-1 , j+1] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i-1 , j+1,layer]
      grid[i-1 , j+1,layer] <- hold
    }
    if (disp == NW) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i-1 , j-1]
      responseTimeGrid[i-1 , j-1] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i-1 , j-1,layer]
      grid[i-1 , j-1,layer] <- hold
    }
    if (disp == SE) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i+1 , j+1]
      responseTimeGrid[i+1 , j+1] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i+1 , j+1,layer]
      grid[i+1 , j+1,layer] <- hold
    }
    if (disp == SW) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i+1 , j-1]
      responseTimeGrid[i+1 , j-1] <<- hold
      hold = grid[i,j,layer]
      grid[i , j,layer] <- grid[i+1 , j-1,layer]
      grid[i+1 , j-1,layer] <- hold
    }
  }
}

spread = function(site, N, NE, NW, E, S, W, SW, SE, probInfect, probReplace, i, j) {
  # SPREAD - Function to return the value of a site
  # at the next time step
  
  if (site == HEALTHY){
    #if site is healthy and has one A1 neighbor, cell becomes infected
    if(N == INFECTEDA1 || NE == INFECTEDA1 || NW == INFECTEDA1 || E == INFECTEDA1 || S == INFECTEDA1 || W == INFECTEDA1 || SW == INFECTEDA1 || SE == INFECTEDA1)
    {
      newSite = INFECTEDA1
      responseTimeGrid[i, j] <<- responseTime
    }
    else if(N+NE+NW+E+W+S+SW+SE > numberOfA2 * 2){
      newSite = INFECTEDA1
    }
    else if(runif(1) < probInfect){
      newSite = INFECTEDA1
    }
    else{
      newSite = HEALTHY
    }
  }
  else if (site == INFECTEDA1){
    #check to see if responseTime has passed for this cell
    if(responseTimeGrid[i, j] == 0){
      newSite = INFECTEDA2
    }else{
      responseTimeGrid[i, j] <<- responseTimeGrid[i, j] - 1
      newSite = INFECTEDA1
    }
  }
  else if (site == INFECTEDA2) {
    newSite = DEAD
  }
  else if (site == DEAD) {
    if(runif(1) < probReplace){
      newSite = HEALTHY
    }else{
      newSite = DEAD
    }
  }
  else if (site == EMPTY) {
    newSite = EMPTY
  }
  return(newSite)
}


### TESTING ###

## test grids = HIV(n, probHIV, probEmpty, probInfect, probReplace, t)
grids = HIV(100, 0.05, 0.4, 0.000001, 0.99, 20)
showGraphs(grids, 15)