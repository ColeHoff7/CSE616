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
  responseTimeGrid <<- periodicLat(responseTimeGrid)
  for (i in 2:(t+1)) {

    bodyExtended = periodicLat(body)
    
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
  newmat1 = latExt
  for (j in 2:(n + 1)) {
    for (i in 2:(n + 1)) {
      site = latExt[i, j]
      if(i > 2){
        N = newmat1[i - 1, j]
      }else{
        N= HEALTHY
      }
      if(i > 2 && j < n){
        NE = newmat1[i-1, j+1]
      }else{
        NE = HEALTHY
      }
      if(i > 2 && j > 2){
        NW = newmat1[i-1, j-1]
      }else{
        NW = HEALTHY
      }
      if(j < n){
        E = newmat1[i, j + 1]
      }else{
        E = HEALTHY
      }
      if(i < n){
        S = newmat1[i + 1, j]
      }else{
        S = HEALTHY
      }
      if(j > 2){
        W = newmat1[i, j - 1]
      }else{
        W = HEALTHY
      }
      if(i < n && j > 2){
        SW = newmat1[i+1, j-1]
      }else{
        SW = HEALTHY
      }
      if(i < n && j < n){
        SE = newmat1[i+1, j+1]
      }else{
        SE = HEALTHY
      }
      if (layer>1){
        
        #if(i < n && i > 2 && j < n && j > 2){
        newmat1= move(site, N, NE, NW, E, S, W, SW, SE, grids, i, j, layer, newmat1, n)
        
        #}
      }
      
    }
  }
  
  newmat2 = matrix(c(rep(0,n*n)), nrow = n)
  for (j in 2:(n + 1)) {
    for (i in 2:(n + 1)) {
      site = newmat1[i, j]
      N = newmat1[i - 1, j]
      NE = newmat1[i-1, j+1]
      NW = newmat1[i-1, j-1]
      E = newmat1[i, j + 1]
      S = newmat1[i + 1, j]
      W = newmat1[i, j - 1]
      SW = newmat1[i+1, j-1]
      SE = newmat1[i+1, j+1]
      newmat2[i - 1, j - 1] = site#spread(site, N, NE, NW, E, S, W, SW, SE, probInfect, probReplace, i, j)
    }
  }
  return(newmat2)
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
    print(length(empty[[1]]))
    plot(healthy[[1]],healthy[[2]],pch=19,col="green",
         xlim=c(0,n+1),ylim=c(0,n+1))
    points(empty[[1]],empty[[2]],col="black",pch=23,bg="white")
    points(infecteda1[[1]],infecteda1[[2]],col="orange",pch=20,bg="orange")
    points(infecteda2[[1]],infecteda2[[2]],col="purple",pch=20,bg="purple")
    points(dead[[1]],dead[[2]],col="black",pch=20,bg="black")
    Sys.sleep(0.2)
  }
}

move = function(site, N, NE, NW, E, S, W, SW, SE, grid, i, j, layer, mat, n) {
  if(site == HEALTHY || site == INFECTEDA1 || site == INFECTEDA2) {
    space = c(0)
    # Finds which spaces around the site are EMPTY and stores them in a vector
    if(N == EMPTY && i>2) {space = c(space, 1)}
    if(NE == EMPTY && i>2 && j<n) {space = c(space, 2)}
    if(NW == EMPTY && i>2 && j>2) {space = c(space, 3)}
    if(E == EMPTY&& j < n) {space = c(space, 4)}
    if(S == EMPTY&&i < n) {space = c(space, 5)}
    if(W == EMPTY&&j>2) {space = c(space, 6)}
    if(SE == EMPTY&&i<n && j>2) {space = c(space, 7)}
    if(SW == EMPTY&&i<n&&j<n) {space = c(space, 8)}
    disp = sample(space,1)
    
    if (disp == 1) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i-1 , j]
      responseTimeGrid[i-1 , j] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i-1 , j]
      mat[i-1 , j] <- hold
    }
    if (disp == 2) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i-1, j+1]
      responseTimeGrid[i-1, j+1] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i-1, j+1]
      mat[i-1, j+1] <- hold
    }
    if (disp == 3) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i-1, j-1]
      responseTimeGrid[i-1, j-1] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i-1, j-1]
      mat[i-1, j-1] <- hold
    }
    if (disp == 4) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i, j + 1]
      responseTimeGrid[i, j + 1] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i, j + 1]
      mat[i, j + 1] <- hold
    }
    if (disp == 5) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i + 1, j]
      responseTimeGrid[i + 1, j] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i + 1, j]
      mat[i + 1, j] <- hold
    }
    if (disp == 6) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i, j - 1]
      responseTimeGrid[i, j - 1] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i, j - 1]
      mat[i, j - 1] <- hold
    }
    if (disp == 7) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i+1, j-1]
      responseTimeGrid[i+1, j-1] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i+1, j-1]
      mat[i+1, j-1] <- hold
    }
    if (disp == 8) {
      hold = responseTimeGrid[i,j]
      responseTimeGrid[i , j] <<- responseTimeGrid[i+1, j+1]
      responseTimeGrid[i+1, j+1] <<- hold
      hold = mat[i,j]
      mat[i , j] <- mat[i+1, j+1]
      mat[i+1, j+1] <- hold
    }
  }
  return(mat)
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
grids = HIV(100, 0.05, .2, 0.000001, 0.99, 10)
showGraphs(grids, 100)