surface = function(expr, from.x=0, to.x=1, from.y=0, to.y=1, n.x=30, 
  n.y=30, col.list=rainbow(30), theta=5, phi=25, mar=c(1,1,1,1), ...) {
  # Build the 2d grid
  x = seq(from=from.x,to=to.x,length.out=n.x)
  y = seq(from=from.y,to=to.y,length.out=n.y)
  plot.grid = expand.grid(x=x,y=y)
  
  # Evaluate the expression to get matrix of z values
  uneval.expr = substitute(expr)
  z.vals = eval(uneval.expr,envir=plot.grid)
  z = matrix(z.vals,nrow=n.x)
  
  # Figure out margins
  orig.mar = par()$mar # Save the original margins
  par(mar=mar)
  col.grid = col.list[round((z[-1,-1]-min(z))/(max(z)-min(z))
    * (length(col.list)-1) + 1)]
  r = persp(x,y,z,theta=theta,phi=phi,col=col.grid,...)
  par(mar=orig.mar) # Restore the original margins
  invisible(r) # Return the persp object invisibly
}
