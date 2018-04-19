HEALTHY = 0
INFECTEDA1 = 1
INFECTEDA2 = 2
DEAD = 3
numberOfA2 = 5
responseTime  = 4 #4 weeks
responseTimeGrid <- array(numeric(), c(0,0))
respoProbTime = 0
therapyTime = 300

HIV = function(n, probHIV, probInfect, probReplace, rank, t) {
  #rank is number from 0 to 8 for effectiveness of drug treatment
  body= initBody( n, probHIV)
  responseTimeGrid <<- body * responseTime
  #Uncomment for step respond prob. for project 2 part b
  #respoProbTime = makeRespoProbTime(.2, .8, t)
  #Uncomment for linear respond prob. for project 2 part c
  respoProbTime = rev(seq(0,t,1))/t
  grids = array(dim=c(n,n,t+1))

  
  grids[,,1] = body
  responseTimeGrid <<- periodicLat(responseTimeGrid)
  prog = .1
  for (i in 2:(t+1)) {
    if(prog*(t+1) < i){
      cat(prog*100,"%  ")
      prog = prog + .1
    }
    bodyExtended = periodicLat(body)
    if(t > therapyTime){
      mode = 1
    }else{
      mode = 0
    }
    body = applyExtended(bodyExtended, probInfect, probReplace, respoProbTime[i], mode, rank)
    grids[,,i] = body
  }
  
  return(grids)
}

makeRespoProbTime = function(min, max, n){
  
  respoProbTime = seq(0,n,1)
  respoProbTime[1:n/2] = max
  respoProbTime[(n/2):n] = min
  return(respoProbTime)
}

applyExtended = function(latExt, probInfect, probReplace, respoProbTime, mode, rank) {
  
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
      newmat[i-1, j-1] = spread(site, N, NE, NW, E, S, W, SW, SE, probInfect, probReplace, i, j, mode, respoProbTime, rank)
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
  savePath = file.choose()
  savePath = substr(savePath,1,nchar(savePath)-4)
  m = dim(graphList)[3]
  for (k in 1:m) {
    #uncomment to save image
    jpeg(paste(savePath,paste(toString(k),".jpg",sep = ""),sep = ""))
    g = graphList[,,k]
    healthy = pointsForGrid(g, HEALTHY)
    infecteda1 = pointsForGrid(g,INFECTEDA1)
    infecteda2 = pointsForGrid(g,INFECTEDA2)
    dead = pointsForGrid(g,DEAD)
    pointSize = 1
    pointType = '.'
    plot(healthy[[1]],healthy[[2]],pch=pointType,col="green",
         xlim=c(0,n+1),ylim=c(0,n+1), cex = pointSize)
    points(infecteda1[[1]],infecteda1[[2]],col="orange",pch=pointType,bg="orange", cex = pointSize)
    points(infecteda2[[1]],infecteda2[[2]],col="purple",pch=pointType,bg="purple", cex = pointSize)
    points(dead[[1]],dead[[2]],col="black",pch=pointType,bg="black", cex = pointSize)
    Sys.sleep(0.2)
  	#paired with jpeg
    dev.off()
  }
}
healthyNum = vector(length = 100)
showTimeGraph = function(graphList,n){
  #Function used to extract number of health, A1, A2, and
  #dead cell and plot vs weeks
  m = dim(graphList)[3]
 # m = dim(grids)[3]
  healthyNum = vector(length = n)
  a1Num = vector(length = n)
  a2Num = vector(length = n)
  deadNum = vector(length = n)
  for (k in 1:m) {
   g = graphList[,,k]
   #For debug
   #g = grids[,,k]
   healthy = pointsForGrid(g, HEALTHY)
   infecteda1 = pointsForGrid(g,INFECTEDA1)
   infecteda2 = pointsForGrid(g,INFECTEDA2)
   dead = pointsForGrid(g,DEAD)
   healthyNum[k] = length(healthy[[1]])
   a1Num[k] = length(infecteda1[[1]])
   a2Num[k] = length(infecteda2[[1]])
   deadNum[k] = length(dead[[1]])
  }
  weeks = seq(1,m,1)
  plot(weeks,healthyNum,type="n",xlab="Weeks",ylab="Num",main = "Cell Infection in Time",col=c("red"))
  legend("topright", legend=c("Healthy","A1Infected","A2Infected","Dead"), lty=c("solid","solid","solid","solid","solid"),col=c("green","orange","red","black"))
  lines(weeks,healthyNum,lwd=2,lty="solid",col=c("green"))
  lines(weeks,a1Num,lwd=2,lty="solid",col=c("orange"))
  lines(weeks,a2Num,lwd=2,lty="solid",col=c("red"))
  lines(weeks,deadNum,lwd=2,lty="solid",col=c("black"))
  return(list(weeks,healthyNum,a1Num,a2Num,deadNum))
}

FFT = function(TimeObject, fs, fRange){
  weeks = TimeObject[[1]]
  healthyNum = TimeObject[[1]]
  a1Num = TimeObject[[1]]
  a2Num = TimeObject[[1]]
  deadNum = TimeObject[[1]]
  
  healthyFFT = fft(healthyNum)
  a1FFT = fft(a1Num)
  a2FFT = fft(a2Num)
  deadFFT = fft(deadNum)
  
  plot.frequency.spectrum(healthyFFT, xlimits=c(0,fRange), fs, "Healthy", "green")
  plot.frequency.spectrum(a1FFT, xlimits=c(0,fRange), fs, "A1", "orange")
  plot.frequency.spectrum(a2FFT, xlimits=c(0,fRange), fs, "A2", "red")
  plot.frequency.spectrum(deadFFT, xlimits=c(0,fRange), fs, "Dead", "black")
  
  
}

#Orignal Author Jo?o Neto++++++++++++++++++++++++++++++++++++++++++++++++++++
#Modified by Wilson Guo++++++++++++++++++++++++++++++++++++++++++++++++++++++
plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k)), fs, label, color) {
  
  zeroPad = seq(0,100,1)
  zeroPad[1:101] = 0
  X.k = fft(c(zeroPad,healthyNum,zeroPad))
  freq = seq(0,(length(X.k)/2),1)*fs/length(X.k)#Construct Freqeuency Axis
  data = abs(X.k/length(test))#Construct Amplitude Axis
  data[2:length(X.k)] = 2*data[2:length(X.k)]#Scale Amplitude
  data = data[1:(length(data)/2+1)]#Cut negative freq
  plot.data  = cbind(freq, data)#Combine Freqeuency and Amplitude
  plot(plot.data, t="h", lwd=2, main="", 
       xlab="Frequency (Hz)", ylab="Strength", 
       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))),col=c(color))
  legend("topright", legend=c(label), lty=c("solid"),col=c(color))
  
}

# Plot the i-th harmonic
# Xk: the frequencies computed by the FFt
#  i: which harmonic
# ts: the sampling time points
# acq.freq: the acquisition rate
plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
  Xk.h <- rep(0,length(Xk))
  Xk.h[i+1] <- Xk[i+1] # i-th harmonic
  harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
  points(ts, harmonic.trajectory, type="l", col=color)
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

spread = function(site, N, NE, NW, E, S, W, SW, SE, probInfect, probReplace, i, j, mode, respoProb, rank) {
  # SPREAD - Function to return the value of a site
  # at the next time step
  # mode : enable probalistic response - 0 for default, 1 to enable
  # respoProb : between 0 and 1
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
    if(mode == 0){
      #check to see if responseTime has passed for this cell
      if(responseTimeGrid[i, j] == 0){
        newSite = INFECTEDA2
      }else{
        responseTimeGrid[i, j] <<- responseTimeGrid[i, j] - 1
        newSite = INFECTEDA1
      }
    }
    else{
      #probalistic promotion from A1 to A2
      if(runif(1) < (1-respoProb)*(rank/8)){
        newSite = INFECTEDA2
      }
      else{
        newSite = INFECTEDA1
      }
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

## test grids = HIV(n, probHIV, probInfect, probReplace, rank, t)
time = 6
num = 400
grids = HIV(num, .05, 0.00001, 0.99, 6, time)
showGraphs(grids, num)
timeCell = showTimeGraph(grids, time+1)
#FFT(timeCell,1,.5)
