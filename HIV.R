HEALTHY = 0
INFECTEDA1 = 1
INFECTEDA2 = 2
DEAD = 3
numberOfA2 = 5
responseTime  = 4*7 #4 weeks
responseTimeGrid <- array(numeric(), c(0,0))

HIV = function(n, probHIV, probInfect, probReplace, t) {
  # FIRE simulation
  body= initBody( n, probHIV)
  responseTimeGrid <<- body * responseTime
  
  grids = array(dim=c(n,n,t+1))

  
  grids[,,1] = body
  responseTimeGrid <<- periodicLat(responseTimeGrid)
  for (i in 2:(t+1)) {
    print(responseTimeGrid)
    bodyExtended = periodicLat(body)
    
    body = applyExtended(bodyExtended, probInfect, probReplace)
    grids[,,i] = body
  }
  
  return(grids)
}

applyExtended = function(latExt, probInfect, probReplace) {
  
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
      newmat[i-1, j-1] = spread(site, N, NE, NW, E, S, W, SW, SE, probInfect, probReplace, i, j)
    }
  }
  return(newmat)
}

initBody = function( n, probHIV ){
  #INITBODY returns an n-by-n grid of values?\
  # HEALTHY (healthy cell), 
  # INFECTEDA1 (infected cell), or INFECTEDA2 (infected cell), or DEAD (dead cell)
  # Pre:	n is the size (number of rows or columns) of the square grid and is positive.
  # 	probHIV is the probability that a cell is initially infected 
  # Post:	 A grid as described above was returned.

  # 1 where infected
  infected = matrix(runif(n^2) < probHIV,nrow=n,ncol=n)
  # 1 where healthy
  healthy = 1 - infected

  body =  healthy * HEALTHY + infected * INFECTEDA1
  
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
    plot(healthy[[1]],healthy[[2]],pch=19,col="green",
         xlim=c(0,n+1),ylim=c(0,n+1))
    points(infecteda1[[1]],infecteda1[[2]],col="red",pch=23,bg="orange")
    points(infecteda2[[1]],infecteda2[[2]],col="red",pch=23,bg="purple")
    points(dead[[1]],dead[[2]],col="red",pch=23,bg="black")
    Sys.sleep(0.2)
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
      responseTimeGrid[i, j] <<- responseTime
    }
    else if(runif(1) < probInfect){
      newSite = INFECTEDA1
      responseTimeGrid[i, j] <<- responseTime
    }
    else{
      newSite = HEALTHY
      responseTimeGrid[i, j] <<- 0
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
    responseTimeGrid[i, j] <<- 0
  }
  else if (site == DEAD) {
    if(runif(1) < probReplace){
      newSite = HEALTHY
      responseTimeGrid[i, j] <<- 0
    }else{
      newSite = DEAD
      responseTimeGrid[i, j] <<- 0
    }
  }
  return(newSite)
}


### TESTING ###

## test grids = HIV(n, probHIV, probInfect, probReplace, t)
grids = HIV(20, .05, 0.00001, 0.99, 50)
showGraphs(grids, 20)