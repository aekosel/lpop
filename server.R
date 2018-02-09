######################################################################
###Author: aekosel
###Purpose: This is the server code for R Shiny. It also requires
###the UI code for use.
######################################################################

#############################
###load libraries and data###
#############################
require(shiny)
require(KernSmooth)
require(nlme)
dat = read.csv("simulated_data.csv")


###################################
###code to find the neighborhood###
###################################
make.box <- function(my.point, n, X, Y, beta.hat=rep(1, 2)){
	#remove NAs
	all.cc = complete.cases(X, Y)
	X = X[all.cc==T, ]
	Y = Y[all.cc==T]

	#marginal scaling
	beta.hat = abs(beta.hat)
	data.mat = sweep(X, 2, beta.hat, '*')
	sc.point = my.point * beta.hat
  
	#find the sup norm distance
	d.supnorm = cbind( abs(data.mat[ , 1] - sc.point[1]), 
		abs(data.mat[ , 2] - sc.point[2]) )
	d.supnorm = apply( d.supnorm, 1, max)

	#choose the n closest points
	dSort = sort( d.supnorm )
	max.d = dSort[n]

	keep = as.logical( d.supnorm <= (max.d + 1e-14) )
	m = sum(keep)
	if( m==n ){
		in.box.supnorm = sort( order(d.supnorm)[1:n] )
	} else {
		keep = d.supnorm <= (max.d + 1e-14)
		m = sum(keep)
		in.box.supnorm = sort( order(d.supnorm)[1:m] )
		d.l1 = abs(data.mat[ , 1] - sc.point[1]) + 
			abs(data.mat[ , 2] - sc.point[2])
		d.break = sort( d.l1[ in.box.supnorm ] )[n] + 1e-14
		keep = (d.supnorm <= (max.d + 1e-14) ) & ( d.l1 <= d.break )
		index = c( 1:length(d.supnorm) )
		in.box.supnorm = index[ keep ]
	}

	#find the corners of the box
	box.corners.supnorm = c( min(X[in.box.supnorm, 1]), max(X[in.box.supnorm, 1]),
    min(X[in.box.supnorm, 2]), max(X[in.box.supnorm, 2]) )
 	names(box.corners.supnorm) = c("x.min", "x.max", "y.min", "y.max")

	#find the mean distance of the box
	d.avg.supnorm = mean( sqrt ( apply( (my.point -  X[in.box.supnorm, ])^2 , 1, sum) ) )
	sc.d.avg.supnorm = mean( sqrt ( apply( (sc.point -  data.mat[in.box.supnorm, ])^2 , 1, sum) ) )
	names(d.avg.supnorm) = "d.avg.supnorm"
	names(sc.d.avg.supnorm) = "sc.d.avg.supnorm"

	#find the balance of the box
	bal.val.supnorm = apply( - my.point +  X[in.box.supnorm, ] , 2, mean)
	sc.bal.val.supnorm = apply( - sc.point +  data.mat[in.box.supnorm, ] , 2, mean)
	bal.no.supnorm = apply( sign( - my.point +  X[in.box.supnorm, ]), 2, mean)

	return( list(
    "in.box"=in.box.supnorm,                          #[[1]] = points in box
    "box.corners"=box.corners.supnorm,                #[[2]] = corners of box
    "d.avg"=c(d.avg.supnorm, sc.d.avg.supnorm),       #[[4]] = mean d of box
    "bal.val"=c(bal.val.supnorm, sc.bal.val.supnorm), #[[5]] = balance of box
    "bal.no"=bal.no.supnorm,                          #[[6]] = balance of box
    "beta.hat"=beta.hat                               #[[7]] = beta hat
     ) )
}


##########################
###code to plot the box###
##########################
plot.box <- function(my.point, n, X, Y, beta.hat=rep(1,2),
	plot.supnorm = T, x1.l=c(0,1), x2.l=c(0,1)){
	#obtain the neighborhood
	out = make.box(my.point, n, X, Y, beta.hat)

	#remove NAs
	all.cc = complete.cases(X, Y)
	X = X[all.cc==T, ]
	Y = Y[all.cc==T]

	#display the point of interest
	plot(my.point[1], my.point[2],
		col="#FF0000", pch=19,
		xlim=x1.l, ylim=x2.l,
		xlab="X1", ylab="X2")

	#display the data
	points(X[ , 1], X[ , 2], pch=19)

	#display the neighborhood
	if(plot.supnorm == T){
		in.box.supnorm = out$in.box
		box.corners.supnorm = out$box.corners
		polygon(c(box.corners.supnorm[1:2], box.corners.supnorm[2:1]),
			c(rep(box.corners.supnorm[3], 2), rep(box.corners.supnorm[4], 2)),
			lty=1, lwd=2, border="#77D3F1")
		points(X[in.box.supnorm, 1], X[in.box.supnorm, 2],
			col="#77D3F1", pch=20)
	}

	#ensure the point of interest is on top
	points(my.point[1], my.point[2],
		col="#FF0000", pch=19)

  return(out)
}


################
###shiny code###
################
shinyServer(function(input, output) {

	###############################################
	###set plot limits based on variables chosen###
	###############################################
	output$x1.lim <- renderUI({
		x1.values <- dat[ , which(names(dat)==input$myvar1)]
		x1.values = x1.values[!is.na(x1.values)]
		sliderInput("x1.lim", "X1 Plot Limits",
			min = floor(min(x1.values)), max = ceiling(max(x1.values)),
			value = c(min(x1.values), max(x1.values)))
	})
	
	output$x2.lim <- renderUI({
		x2.values <- dat[ , which(names(dat)==input$myvar2)]
		x2.values = x2.values[!is.na(x2.values)]
		sliderInput("x2.lim", "X2 Plot Limits",
			min = floor(min(x2.values)), max = ceiling(max(x2.values)),
			value = c(min(x2.values), max(x2.values)))
	})

	########################################################
	###set up slider for number of points in neighborhood###
	########################################################
	output$n <- renderUI({
	all.cc = complete.cases(dat[ , which(names(dat)==input$myout)],
		dat[ , which(names(dat)==input$myvar1)],
		dat[ , which(names(dat)==input$myvar2)])
	nopt = length(all.cc[all.cc==T])
	sliderInput("n", "No. of Points in Box:", 2, nopt, 0.1*nopt, step = 1,
		animate=animationOptions(interval=100, loop=F))  
	})
     
	##############################
	###display the neighborhood###
	##############################
	output$myplot <- renderPlot({
		Yv = dat[ , which(names(dat)==input$myout)]
		Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
		all.cc = complete.cases(Xv, Yv)
		Xv = Xv[all.cc==T, ]
		Yv = Yv[all.cc==T]

	if(input$scaled == "No"){
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, beta.hat=rep(1,2))
	}
	else if(input$scaled == "Globally"){
		my.beta.hat = rep(1,2)
		my.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
		my.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
    }
	else{
		my.beta.hat = rep(1,2)
		my.band.1 = dpill(Xv[, 1], Yv)
		beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
		my.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
		my.band.2 = dpill(Xv[, 2], Yv)
		beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
		my.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
    
    #set up the plot window
	plot(input$my.point.x1, input$my.point.x2,
		col="#FF0000", pch=19,
		xlim=input$x1.lim, ylim=input$x2.lim,
		xlab=input$myvar1, ylab=input$myvar2,
		main="Covariate Distribution")

	#display the data
	points(Xv[ , 1], Xv[ , 2],
	pch=19)

    #display the neighborhood
	in.box.supnorm = my.box.vals$in.box
	box.corners.supnorm = my.box.vals$box.corners
	polygon(c(box.corners.supnorm[1:2], box.corners.supnorm[2:1]),
	c(rep(box.corners.supnorm[3], 2), rep(box.corners.supnorm[4], 2)),
		lty=1, lwd=2, border="#77D3F1")
	points(Xv[in.box.supnorm, 1], Xv[in.box.supnorm, 2],
		col="#77D3F1", pch=20)

	#display the point of interest
	points(input$my.point.x1, input$my.point.x2,
		col="#FF0000", pch=19)
	})

	################################
	###make the outcome histogram###
	################################
	output$myhist <- renderPlot({
		Yv = dat[ , which(names(dat)==input$myout)]
		Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
		all.cc = complete.cases(Xv, Yv)
		Xv = Xv[all.cc==T, ]
		Yv = Yv[all.cc==T]

	if(input$scaled == "No"){
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, beta.hat=rep(1,2))
	}
	else if(input$scaled == "Globally"){
		my.beta.hat = rep(1,2)
		my.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
		my.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	else{
		my.beta.hat = rep(1,2)
		my.band.1 = dpill(Xv[, 1], Yv)
		beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
		my.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
		my.band.2 = dpill(Xv[, 2], Yv)
		beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
		my.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	ib.m = quantile(Yv[my.box.vals$in.box], input$histquant/100)

	#plot the histogram
	hist(Yv[my.box.vals$in.box], col="#77D3F1",
		xlab="Outcomes", main="Histogram of Neighbors",
		xlim=c(min(Yv),max(Yv)),
		breaks=input$histboxes)
	text((min(Yv) + max(Yv))/2,
		max(hist(Yv[my.box.vals$in.box], col="#77D3F1",
		xlab="Outcomes", main="Histogram of Neighbors",
		xlim=c(min(Yv),max(Yv)),
		breaks=input$histboxes)$counts),
		paste("Percentile", input$histquant, "is", ib.m, sep=" "))
		abline(v=ib.m, col="red", lwd=2)
	})

	####################################
	###calculate box scale parameters###
	####################################
	scaleValues <- reactive({
	#obtain the neighborhood values
	Yv = dat[ , which(names(dat)==input$myout)]
	Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
	all.cc = complete.cases(Xv, Yv)
	Xv = Xv[all.cc==T, ]
	Yv = Yv[all.cc==T]

	global.beta.hat = rep(1,2)
	global.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])),
		crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
	global.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])),
		crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]

	local.beta.hat = rep(1,2)
	my.band.1 = dpill(Xv[, 1], Yv)
	beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
	local.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
	my.band.2 = dpill(Xv[, 2], Yv)
	beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
	local.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
      
	my.data = data.frame(c(1, global.beta.hat[1], local.beta.hat[1]),
		c(1, global.beta.hat[2], local.beta.hat[2]),
		row.names=c("No", "Globally", "Locally"))
	names(my.data) = c(paste(input$myvar1," Scale", sep=""), paste(input$myvar2," Scale", sep=""))
	my.data = my.data
	})

	####################################
	###calculate confidence intervals###
	####################################
	ciValues <- reactive({
	Yv = dat[ , which(names(dat)==input$myout)]
	Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
	all.cc = complete.cases(Xv, Yv)
	Xv = Xv[all.cc==T, ]
	Yv = Yv[all.cc==T]

	#calculate means and errors
	if(input$scaled == "No"){
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, beta.hat=rep(1,2))
	}
	else if(input$scaled == "Globally"){
		my.beta.hat = rep(1,2)
		my.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
		my.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	else{
		my.beta.hat = rep(1,2)
		my.band.1 = dpill(Xv[, 1], Yv)
		beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
		my.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
		my.band.2 = dpill(Xv[, 2], Yv)
		beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
		my.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	ib.m = mean(Yv[my.box.vals$in.box])
	ib.se = sd(Yv[my.box.vals$in.box])/input$n
	o.m = mean(Yv)
	o.se = sd(Yv)/length(Yv)
   
	r.1 = c(ib.m - 1.96*ib.se, ib.m + 1.96*ib.se)
	r.2 = c(o.m - 1.96*o.se, o.m + 1.96*o.se) 
	my.data = data.frame(rbind(r.1, r.2), row.names=c("in box ci", "overall ci"))
	names(my.data) = c("lower", "upper")  
	my.data = my.data
	})

	#########################
	###calculate quantiles###
	#########################
	quantileValues <- reactive({
	Yv = dat[ , which(names(dat)==input$myout)]
	Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
	all.cc = complete.cases(Xv, Yv)
	Xv = Xv[all.cc==T, ]
	Yv = Yv[all.cc==T]

	#calculate means and errors
	if(input$scaled == "No"){
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, beta.hat=rep(1,2))
	}
	else if(input$scaled == "Globally"){
		my.beta.hat = rep(1,2)
		my.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
		my.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	else{
		my.beta.hat = rep(1,2)
		my.band.1 = dpill(Xv[, 1], Yv)
		beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
		my.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
		my.band.2 = dpill(Xv[, 2], Yv)
		beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
		my.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	ib.quant = quantile(Yv[my.box.vals$in.box], c(.25, .5, .75))
	o.quant = quantile(Yv, c(.25, .5, .75))

	my.data = data.frame(cbind(ib.quant, o.quant), row.names=c("25th", "50th", "75th"))
	names(my.data) = c("in box", "overall")  
	my.data = my.data
	})

	####################################
	###calculate neighborhood corners###
	####################################
	cornerValues <- reactive({
	Yv = dat[ , which(names(dat)==input$myout)]
	Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
	all.cc = complete.cases(Xv, Yv)
	Xv = Xv[all.cc==T, ]
	Yv = Yv[all.cc==T]
	if(input$scaled == "No"){
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, beta.hat=rep(1,2))
	}
	else if(input$scaled == "Globally"){
		my.beta.hat = rep(1,2)
		my.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
		my.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])),
			crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	else{
		my.beta.hat = rep(1,2)
		my.band.1 = dpill(Xv[, 1], Yv)
		beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
		my.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
		my.band.2 = dpill(Xv[, 2], Yv)
		beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
		my.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	my.data = data.frame( cbind(c(my.box.vals$box.corners[1],my.box.vals$box.corners[3]),
		c(my.box.vals$box.corners[2],my.box.vals$box.corners[4])),
		row.names=c(input$myvar1, input$myvar2) )
		names(my.data) = c("min","max")

	my.data = my.data
	})

	##########################################
	###count the points in the neighborhood###
	##########################################
	inboxValues <- reactive({
	#get the box values
	Yv = dat[ , which(names(dat)==input$myout)]
	Xv = cbind(dat[ , which(names(dat)==input$myvar1)], dat[ , which(names(dat)==input$myvar2)])
	all.cc = complete.cases(Xv, Yv)
	Xv = Xv[all.cc==T, ]
	Yv = Yv[all.cc==T]
	if(input$scaled == "No"){
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, beta.hat=rep(1,2))
	}
	else if(input$scaled == "Globally"){
		my.beta.hat = rep(1,2)
		my.beta.hat[1] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 1])), crossprod(cbind(rep(1,length(Yv)), Xv[ , 1]), Yv) )[2]
		my.beta.hat[2] = solve(crossprod(cbind(rep(1,length(Yv)), Xv[ , 2])), crossprod(cbind(rep(1,length(Yv)), Xv[ , 2]), Yv) )[2]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	else{
		my.beta.hat = rep(1,2)
		my.band.1 = dpill(Xv[, 1], Yv)
		beta.est1 = locpoly(Xv[, 1], Yv, drv=1, degree=2, bandwidth=my.band.1)
		my.beta.hat[1] = beta.est1$y[which(abs(beta.est1$x-input$my.point.x1) == min(abs(beta.est1$x-input$my.point.x1)))]
 		my.band.2 = dpill(Xv[, 2], Yv)
		beta.est2 = locpoly(Xv[, 2], Yv, drv=1, degree=2, bandwidth=my.band.2)
		my.beta.hat[2] = beta.est2$y[which(abs(beta.est2$x-input$my.point.x2) == min(abs(beta.est2$x-input$my.point.x2)))]
		my.box.vals = make.box(c(input$my.point.x1, input$my.point.x2), input$n, Xv, Yv, my.beta.hat)
	}
	pts.keep = length(my.box.vals$in.box)

	#put the counts in a table
	my.data = data.frame("Actual"=pts.keep, "Requested"=as.integer(input$n), row.names="Points")
	})


	###########################################
	###put all the values into output tables###
	###########################################
	output$boxsizetable <- renderTable({
		scaleValues()
	})

	output$citable <- renderTable({
		ciValues()
	})

	output$quantiletable <- renderTable({
		quantileValues()
	})

	output$cornertable <- renderTable({
		cornerValues()
	})

	output$inboxtable <- renderTable({
		inboxValues()
	})

})

