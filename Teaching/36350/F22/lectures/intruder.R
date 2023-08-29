intruder.sim = function(n, p, maxsteps=1000, dt=0.5, seed=NULL, plot=TRUE, 
                        verbose=FALSE) {
  
  if (!is.null(seed)) set.seed(seed)
  if (!plot) dt=0
  
  # Intruder's initial location
  x = sample(1:n,1)
  y = sample(1:n,1)
  
  t = 1
  A = NULL
  
  while (t<=maxsteps) {
    # Generate noise 
    X = matrix(as.numeric(runif(n^2)<=p),n,n)
    
    # Intruder step
    d = random.step(x,y,n)
    x = x+d[1]
    y = y+d[2]
    X[x,y] = 1
    
    # Run our algorithm
    A = find.candidates(X,A)
    
    if (plot) {
      plot.onetime(X,t,identify=TRUE,A,x,y)
    }
    
    if (nrow(A)==1) {
      break
    }
    else if (verbose) {
      cat(sprintf("At t=%i, we have %i candidates...\n",t,nrow(A)))
    }
    
    Sys.sleep(dt)
    t = t+1
  }
  
  if (nrow(A)>1) {
    t = t-1
  }
  else if (verbose) {
    cat(sprintf("At t=%i, we caught him! The intruder is at (%i,%i).\n",
                t,A[1,1],A[1,2]))
  }
  
  return(t)
}

random.step = function(x, y, n) {
  if (1<x & x<n & 1<y & y<n) {
    i = sample(1:5,1)
    return(matrix(c(-1,1,0,0,0,0,0,0,-1,1),nrow=5)[i,])
  }
  else if (x==1 & 1<y & y<n) {
    i = sample(1:4,1)
    return(matrix(c(1,0,0,0,0,0,-1,1),nrow=4)[i,])
  }
  else if (x==n & 1<y & y<n) {
    i = sample(1:4,1)
    return(matrix(c(-1,0,0,0,0,0,-1,1),nrow=4)[i,])
  }
  else if (1<x & x<n & y==1) {
    i = sample(1:4,1)
    return(matrix(c(-1,1,0,0,0,0,0,1),nrow=4)[i,])
  }
  else if (1<x & x<n & y==n) {
    i = sample(1:4,1)
    return(matrix(c(-1,1,0,0,0,0,0,-1),nrow=4)[i,])
  }
  else {
    i = sample(1:3,1)
    return(matrix(c(c(-1,1)[(x==1)+1],0,0,0,0,c(-1,1)[(y==1)+1]),
                  nrow=3)[i,])
  }
}

plot.onetime = function(X, t, identify=TRUE, A, x0, y0, ...) {
  n=nrow(X)
  par(mar=c(1,1,3,1))
  if (n<=30) pch=19
  else pch=20
  main=sprintf("t = %i",t)
  if (identify) {
    main=sprintf("%s\n(number of candidates = %i)",main,nrow(A))
  }
  lim = c(0,n)
  plot(c(),c(),xlim=lim,ylim=lim,axes=FALSE,main=main,xlab="",ylab="",...)
  
  if (identify) {
    rect(A[,1]-1,A[,2]-1,A[,1],A[,2],col="red")
  }
  
  i=which(X!=0)
  x=(i-1)%%n+1
  y=ceiling(i/n)
  points(x-0.5,y-0.5,col="black",pch=pch)
  
  if (identify) {
    points(x0-0.5,y0-0.5,col="white",pch=pch)
  }
  
  segments(rep(0,n+1),0:n,rep(n,n+1),0:n,col="gray")
  segments(0:n,rep(0,n+1),0:n,rep(n,n+1),col="gray")
}

# This is the actual algorithmic part
find.candidates = function(X, A) {
  i = which(X!=0)
  AA = matrix(nrow=length(i),ncol=2)
  n = nrow(X)
  AA[,1] = (i-1)%%n+1
  AA[,2] = floor((i-1)/n)+1
  
  # If we've only seen one timepoint, everything is a candidate
  if (is.null(A)) return(AA)
  
  # Otherwise, refine our list by checking neighboring spots
  ii = rep(FALSE,nrow(AA))
  for (j in 1:nrow(AA)) {
    ii[j] = min(rowSums(abs(scale(A,AA[j,],FALSE)))) <= 1
  }
  return(AA[ii,,drop=FALSE])
}

check.neighbors = function(X, x, y) {
  n = nrow(X)
  if (1<x & x<n & 1<y & y<n) {
    a = matrix(c(x-1,x+1,x,x,x,y,y,y,y-1,y+1),nrow=5)
  }
  else if (x==1 & 1<y & y<n) {
    a = matrix(c(x+1,x,x,x,y,y,y-1,y+1),nrow=4)
  }
  else if (x==n & 1<y & y<n) {
    a = matrix(c(x-1,x,x,x,y,y,y-1,y+1),nrow=4)
  }
  else if (1<x & x<n & y==1) {
    a = matrix(c(x-1,x+1,x,x,y,y,y,y+1),nrow=4)
  }
  else if (1<x & x<n & y==n) {
    a = matrix(c(x-1,x+1,x,x,y,y,y,y-1),nrow=4)
  }
  else {
    a = matrix(c(x+c(-1,1)[(x==1)+1],x,x,y,y,y+c(-1,1)[(y==1)+1]),nrow=3)
  }
  return(a[X[a]!=0,])
}